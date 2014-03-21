library(plyr)
library(MCMCglmm)
library(gdata)
library(matrixcalc)
library(ggplot2)
library(Morphometrics)
library(reshape2)
library(doMC)
library(foreach)

load("./maindata.RData")

getModel <- function(otu){
  models = otu$model[!is.na(otu$model)]
  models = models[[1]]
  fixed = as.character(models$call$formula)[[3]]
  if(fixed == "1") return(NULL)
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
  if(is.null(otu$fixed)){
    current.matrix = otu$vcv[paste("D", num.traits, sep='')][[1]]
    trait.names = rownames(current.matrix)
    if(!is.positive.definite (current.matrix)){
      current.matrix = as.matrix(nearPD(current.matrix)[[1]])
      rownames(current.matrix) <- colnames(current.matrix) <- trait.names 
    }
    n = dim(otu$data[complete.cases(otu$data[,trait.names]),])[1]
    Ps = MonteCarloStat(current.matrix, num.traits+1, 1000, ComparisonFunc=function(x, y) y, cov)
    Ps = aperm(Ps, c(2,3,1))
    dimnames(Ps) = list(trait.names, trait.names)
  }
  else{
    trait.names = rownames(otu$vcv[[paste("D", num.traits, sep='')]])
    traits = paste("cbind(", paste(trait.names, collapse=','), ")", sep = '')
    formula = paste(traits, otu$fixed, sep = "~")
    prior = list(R = list(V = diag(num.traits), n = num.traits+1))
    mcmc.model = MCMCglmm( as.formula(formula),
                           data = otu$data[complete.cases(otu$data[,trait.names]),],
                           prior = prior,
                           rcov = ~us(trait):units,
                           family = rep("gaussian", num.traits),
                           verbose = FALSE)
    Ps = array(mcmc.model$VCV, dim = c(1000, num.traits, num.traits))
    Ps = aperm(Ps, c(2,3,1))
    dimnames(Ps) = list(trait.names, trait.names)
  }
  return(Ps)
}

llply(x, getModel)

for(i in 1:length(x)){
  x[[i]]$data <- cbind(x[[i]]$raw, x[[i]]$categorical)
  x[[i]]$fixed <- getModel(x[[i]])
}
x$Myrmecophaga$fixed <- NULL

generateMCMCArray = function(num.traits){
  mask = !laply(x, function(otu) all(is.na(otu$vcv[[paste("D", num.traits, sep='')]])))
  Ps = laply(x[mask], function(x) runMCMCModel(x, num.traits), .progress='text')
  Ps = aperm(Ps, c(2,3,1,4))
  dimnames(Ps)[[3]] = names(x[mask])
  return(Ps)
}

# Ps = generateMCMCArray(35)
# Ps = Ps[,,,901:1000]
# save(Ps, file = 'xenatraMCMC.samples.Rdata')
load('./xenatraMCMC.samples.Rdata')