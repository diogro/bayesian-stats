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
  return(Ps)
}

models = llply(x, getModel)

for(i in 1:length(x)){
  x[[i]]$data <- cbind(x[[i]]$raw, x[[i]]$categorical)
  x[[i]]$fixed <- getModel(x[[i]])
}