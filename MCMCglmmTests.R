# creating some data

data = read.csv("dados.csv")
data$individual = 

# required libs

library(MCMCglmm)
library(gdata)
library(matrixcalc)
library(ggplot2)
library(Morphometrics)
library(reshape2)

# fitting a simple model

lm.model = lm(as.matrix(data[,1:4])~data[,5])

# fitting MCMCglmm model

mcmc.model = MCMCglmm(cbind(UMERO, ULNA)~trait:ESPECIE - 1 , data = data, family = rep("gaussian", 2))

