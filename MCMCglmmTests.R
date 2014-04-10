library(MCMCglmm)
library(gdata)
library(matrixcalc)
library(ggplot2)
library(Morphometrics)
library(reshape2)

# creating some data

data = read.csv("dados.csv")

num.traits = 4

# fitting a simple model

lm.model = lm(as.matrix(data[,1:num.traits])~data$ESPECIE)

lm.model.Total = lm(as.matrix(data[,1:num.traits])~1)
# fitting MCMCglmm model

prior = list(R = list(V = diag(num.traits), n = 2))

#prior = list(R = list(V = diag(num.traits), n = 2),
             #G = list(G1 = list(V = diag(num.traits)*0.002, n = num.traits+1)))


trait.names = paste("cbind(", paste(names(data)[1:num.traits], collapse=','), ")", sep = '')
#fixed.effects = "trait:as.factor(ESPECIE) - 1"
fixed.effects = "trait - 1"
formula = paste(trait.names, fixed.effects, sep = "~")

mcmc.model = MCMCglmm( as.formula(formula),
                       data = data,
                       prior = prior,
                       rcov = ~us(trait),
                       family = rep("gaussian", num.traits),
                       verbose = TRUE)
G.mcmc = apply(array(mcmc.model$VCV[,1:(num.traits*num.traits)], dim = c(1000, num.traits, num.traits)), 2:3, mean)
Ps = array(mcmc.model$VCV, dim = c(1000, num.traits, num.traits))
avP = apply(Ps, 2:3, mean)
print(avP)
