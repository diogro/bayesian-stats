source('./cov_tensors.R')
source('xenartraSamples.R')


library(ape)
library(Morphometrics)
TREE <- "(((Bradypus:1,Choloepus:1,Paramylodon:1):1,((Tamandua:1,Myrmecophaga:1):1,Cyclopes:1):1):1,(Dasypus:1,(((Tolypeutes:1,Cabassous:1):1,Priodontes:1):1,((Zaedyus:1,Chaetophractus:1):1,Euphractus:1):1):1):1);"
tree.25 <- read.tree(text = TREE)

TREE <- "(((Bradypus:1,Choloepus:1):1,((Tamandua:1,Myrmecophaga:1):1,Cyclopes:1):1):1,(Dasypus:1,(((Tolypeutes:1,Cabassous:1):1,Priodontes:1):1,((Zaedyus:1,Chaetophractus:1):1,Euphractus:1):1):1):1);"
tree.35 <- read.tree(text = TREE)

load('./maindata.RData')
sample.size = ldply(x, function(x) as.data.frame(x$df))

covTensor = function(Ps, tree, sample.size) {
  avgPs = llply(alply(apply(Ps, 1:3, mean), 3), function(x) as.matrix(nearPD(x)[[1]]))
  names(avgPs) = dimnames(Ps)[[3]]
  avgP = AncestralStates(tree, avgPs, sample.size)[[tree$edge[1,1]]]
  
  m = dim(Ps)[3]
  n = dim(Ps)[1]
  rand.Ps = laply(dimnames(Ps)[[3]], 
                  function(x) MonteCarloStat(avgP, n+1, 100, 
                                             ComparisonFunc=function(x, y) y, cov))
  rand.Ps = aperm(rand.Ps, c(3, 4, 1, 2))
  
  MCMC.covtensor <- covtensor(Ps)
  MCMC.covtensor.rand <- covtensor(rand.Ps)
  return(list(obs = MCMC.covtensor,
              rand = MCMC.covtensor.rand))
}

# MCMC.covTensor = list()
# 
# load('./xenartraMCMC.25.samples.Rdata')
# MCMC.covTensor[['25']] = covTensor(Ps, tree.25, sample.size$D25)
# 
# load('./xenartraMCMC.28.samples.Rdata')
# MCMC.covTensor[['28']] = covTensor(Ps, tree.25, sample.size$D28)
# 
# load('./xenartraMCMC.32.samples.Rdata')
# MCMC.covTensor[['32']] = covTensor(Ps, tree.35, sample.size$D32)
# 
# load('./xenartraMCMC.35.samples.Rdata')
# MCMC.covTensor[['35']] = covTensor(Ps, tree.35, sample.size$D35)
# 
# save(MCMC.covTensor, file = 'xenartraCovTensor.Rdata')
load('xenartraCovTensor.Rdata')

dat.cov.tensor = ldply(names(MCMC.covTensor), 
                       function(x) CovTensorDataFrame (MCMC.covTensor[[x]][[1]], 
                                                       MCMC.covTensor[[x]][[2]],
                                                       as.numeric(x)))

cov.tensor.plot = PlotCovTensor(dat.cov.tensor)


load('./xenartraMCMC.35.samples.Rdata')

m = dim(Ps)[3]
n = dim(Ps)[1]
nnonzero <- min(n*(n+1)/2,m-1)

HPD.eT.val <- cbind(HPDinterval(as.mcmc(MCMC.covTensor[[3]][[1]]$MCMCSeigvals[,1:nnonzero]), prob=0.95), 
                    HPDinterval(as.mcmc(MCMC.covTensor[[3]][[2]]$MCMCSeigvals[,1:nnonzero]), prob=0.95))
round(HPD.eT.val,3)
round(MCMC.covTensor[[3]][[1]]$ordered.tensor.summary[1:n*2,2:9],3)


