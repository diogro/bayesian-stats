library(plyr)
library(MCMCglmm)
library(gdata)
library(matrixcalc)
library(ggplot2)
library(Morphometrics)
library(reshape2)

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
    trait.names = rownames(otu$vcv[[paste("D", num.traits, sep='')]])
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
        #prior = list(R = list(V = otu$vcv[[paste("D", num.traits, sep='')]], n = num.traits+50))
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
        vars = mcmcVar(otu, num.traits)
        Ps = Ps * vars
        Ps = aperm(Ps, c(2,3,1))
        dimnames(Ps) = list(trait.names, trait.names)
        return(Ps)
    })
    return(Ps)
}

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
    return(vars)
}

llply(x, getModel)

for(i in 1:length(x)){
  x[[i]]$data <- cbind(x[[i]]$raw, x[[i]]$categorical)
  x[[i]]$fixed <- getModel(x[[i]])
}
x$Myrmecophaga$fixed <- "trait - 1"

generateMCMCArray = function(num.traits){
  mask = !laply(x, function(otu) all(is.na(otu$vcv[[paste0("D", num.traits)]])))
  Ps = laply(x[mask], function(x) runMCMCModel(x, num.traits), .progress='text', .inform = TRUE)
  Ps = aperm(Ps, c(2,3,1,4))
  dimnames(Ps)[[3]] = names(x[mask])
  return(Ps)
}

#Ps = list()
#Ps[['25']] = generateMCMCArray(25)
#Ps[['28']] = generateMCMCArray(28)
#Ps[['32']] = generateMCMCArray(32)
#Ps[['35']] = generateMCMCArray(35)
#save(Ps, file = "xenartraMCMCsamples.Rdata")
load("./xenartraMCMCsamples.Rdata")

#num.traits = 25
#otu = x$Myrmecophaga
#Ps = runMCMCModel(x$Myrmecophaga, 25)
#RandomSkewers(apply(Ps, 1:2, mean), x$Myrmecophaga$vcv$D25)
#Ps = runMCMCModel(x$Myrmecophaga, 28)
#RandomSkewers(apply(Ps, 1:2, mean), x$Myrmecophaga$vcv$D28)
#Ps = runMCMCModel(x$Myrmecophaga, 32)
#RandomSkewers(apply(Ps, 1:2, mean), x$Myrmecophaga$vcv$D32)
#Ps = runMCMCModel(x$Myrmecophaga, 35)
#RandomSkewers(apply(Ps, 1:2, mean), x$Myrmecophaga$vcv$D35)

#avPs = apply(Ps[['25']], 1:3, mean)
#mat.list = alply(avPs, 3)
#mats = llply(x, function(x) x$vcv$D25)
#for(i in 1: length(mats)) { print(names(x)[[i]]); print(RandomSkewers(mats[[i]], mat.list[[i]]))}
