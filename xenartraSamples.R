library(plyr)
library(MCMCglmm)
library(gdata)
library(matrixcalc)
library(ggplot2)
library(Morphometrics)
library(reshape2)
Sys.setlocale(locale="C")

load("./maindata.RData")

getModel <- function(otu){
  models = otu$model[!is.na(otu$model)]
  models = models[[1]]
  fixed = as.character(models$call$formula)[[3]]
  if(fixed == "1") return("trait - 1")
  fixed = gsub("[\\(\\)]", "", regmatches(fixed, gregexpr("\\(.*?\\)", fixed))[[1]])
  if(length(fixed) > 1) {
    fixed.effects = paste(paste(paste("trait:", fixed, sep = ''), collapse= ' * '), " - 1", sep = '')
  }
  else{
    fixed.effects = paste("trait:", fixed, " - 1", sep = '')
  }
  return(fixed.effects)
}

runMCMCModel <- function(otu, num.traits){
    trait.names = rownames(otu$vcv[[paste0("D", num.traits)]])
    traits = paste("cbind(", paste(trait.names, collapse=','), ")", sep = '')
    formula = paste(traits, otu$fixed, sep = "~")
    prior = list(R = list(V = diag(num.traits), n = num.traits+1))
    Ps = tryCatch({
        mcmc.model = MCMCglmm(as.formula(formula),
                              data = otu$data[complete.cases(otu$data[,trait.names]),],
                              prior = prior,
                              rcov = ~us(trait):units,
                              family = rep("gaussian", num.traits),
                              thin = 100,
                              verbose = FALSE)
        Ps = array(mcmc.model$VCV, dim = c(100, num.traits, num.traits))
        Ps = aperm(Ps, c(2,3,1))
        dimnames(Ps) = list(trait.names, trait.names)
        return(Ps)
    }, error = function(cond){
        warning("scaling data")
        print(otu$data$ESPECIE[[1]])
        data = data.frame(scale(otu$data[complete.cases(otu$data[,trait.names]),trait.names], scale = TRUE),
                          otu$data[complete.cases(otu$data[,trait.names]),!names(otu$data) %in% trait.names])
        mcmc.model = MCMCglmm(as.formula(formula),
                              data = data,
                              prior = prior,
                              rcov = ~us(trait):units,
                              thin = 100,
                              family = rep("gaussian", num.traits),
                              verbose = FALSE)
        Ps = array(mcmc.model$VCV, dim = c(100, num.traits, num.traits))
        Ps = aaply(Ps, 1, cov2cor)
        vars = mcmcVar(otu, num.traits)
        #Ps = aaply(Ps, 1, '*', vars)
        Ps = Ps * vars
        Ps = aperm(Ps, c(2,3,1))
        dimnames(Ps) = list(trait.names, trait.names)
        return(Ps)
    })
    return(Ps)
}

otu = x$'Myrmecophaga'
num.traits = 25
mcmcVar <- function(otu, num.traits){
    trait.names = rownames(otu$vcv[[paste("D", num.traits, sep='')]])
    traits = paste("cbind(", paste(trait.names, collapse=','), ")", sep = '')
    formula = paste(traits, otu$fixed, sep = "~")
    prior = list(R = list(V = diag(num.traits), n = num.traits+1))
    data = data.frame(scale(otu$data[complete.cases(otu$data[,trait.names]),trait.names], scale = FALSE),
                      otu$data[complete.cases(otu$data[,trait.names]),!names(otu$data) %in% trait.names])
    mcmc.model = MCMCglmm(as.formula(formula),
                          data = data,
                          prior = prior,
                          rcov = ~idh(trait):units,
                          family = rep("gaussian", num.traits),
                          thin = 100,
                          verbose = FALSE)
    vars = aaply(mcmc.model$VCV, 1, function(x) outer(sqrt(x), sqrt(x)))
    #vars = colMeans(vars)
    return(vars)
}

llply(x, getModel)

for(i in 1:length(x)){
  x[[i]]$data <- cbind(x[[i]]$raw, x[[i]]$categorical)
  x[[i]]$fixed <- getModel(x[[i]])
}


generateMCMCArray = function(num.traits){
  mask = !laply(x, function(otu) all(is.na(otu$vcv[[paste0("D", num.traits)]])))
  Ps = laply(x[mask], function(x) runMCMCModel(x, num.traits), .progress='text', .inform = TRUE, .parallel = TRUE)
  Ps = aperm(Ps, c(2,3,1,4))
  dimnames(Ps)[[3]] = names(x[mask])
  return(Ps)
}

library(doMC)
registerDoMC(10)

Ps = list()
Ps[['25']] = generateMCMCArray(25)
Ps[['28']] = generateMCMCArray(28)
Ps[['32']] = generateMCMCArray(32)
Ps[['35']] = generateMCMCArray(35)
save(Ps, file = "xenartraMCMCsamples.Rdata")
#load("./xenartraMCMCsamples.Rdata")

avPs = apply(Ps[['25']], 1:3, mean)
mat.list = alply(avPs, 3)
mats = llply(x, function(x) x$vcv$D25)
(Map(function(x, y) MatrixCompare(x, y)[,1], mats, mat.list))

i = 1
RandomSkewers(alply(Ps[['25']][,,i,], 3), mats[[i]], num.cores=10)
plot(laply(mats, function(x) sum(diag(x)))~laply(mat.list, function(x) sum(diag(x))))
