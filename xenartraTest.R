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
    Ps = MonteCarloStat(current.matrix, 40, 1000, ComparisonFunc=function(x, y) y, cov)
    Ps = aperm(Ps, c(2,3,1))
    trait.names = rownames(current.matrix)
    dimnames(Ps) = list(trait.names, trait.names)  
  }
  else{
    trait.names = rownames(otu$vcv[paste("D", num.traits, sep='')][[1]])
    trait.names = paste("cbind(", paste(trait.names, collapse=','), ")", sep = '')
    formula = paste(trait.names, otu$fixed, sep = "~")
    prior = list(R = list(V = diag(num.traits), n = 2))
    mcmc.model = MCMCglmm( as.formula(formula),
                           data = otu$data,
                           prior = prior,
                           rcov = ~us(trait):units,
                           family = rep("gaussian", num.traits),
                           verbose = TRUE)
    Ps = array(mcmc.model$VCV, dim = c(1000, num.traits, num.traits))
    Ps = aperm(Ps, c(2,3,1))
  }
  return(Ps)
}

models = llply(x, getModel)

for(i in 1:length(x)){
  x[[i]]$data <- cbind(x[[i]]$raw, x[[i]]$categorical)
  x[[i]]$fixed <- getModel(x[[i]])
}

Ps = laply(x[1:2], function(x) runMCMCModel(x, 35))
Ps = aperm(Ps, c(2,3,1,4))
dimnames(Ps)[[3]] = names(x[1:2])
