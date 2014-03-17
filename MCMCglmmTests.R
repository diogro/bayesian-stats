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

# fitting MCMCglmm model

prior = list(R = list(V = diag(num.traits), n = 2))

mcmc.model = MCMCglmm(cbind(UMERO, ULNA, FEMUR, TIBIA)~trait:ESPECIE - 1 , 
                      data = data, 
                      prior = prior,
                      rcov = ~us(trait):units,
                      family = rep("gaussian", num.traits),
                      verbose = TRUE)

Ps = array(mcmc.model$VCV, dim = c(1000, num.traits, num.traits))
avP = apply(Ps, 2:3, mean)