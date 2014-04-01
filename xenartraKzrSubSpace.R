source('./KrzSubspace.R')

library(Morphometrics)
library(ape)
TREE <- "(((Bradypus:1,Choloepus:1,Paramylodon:1):1,((Tamandua:1,Myrmecophaga:1):1,Cyclopes:1):1):1,(Dasypus:1,(((Tolypeutes:1,Cabassous:1):1,Priodontes:1):1,((Zaedyus:1,Chaetophractus:1):1,Euphractus:1):1):1):1);"
tree.25 <- read.tree(text = TREE)

TREE <- "(((Bradypus:1,Choloepus:1):1,((Tamandua:1,Myrmecophaga:1):1,Cyclopes:1):1):1,(Dasypus:1,(((Tolypeutes:1,Cabassous:1):1,Priodontes:1):1,((Zaedyus:1,Chaetophractus:1):1,Euphractus:1):1):1):1);"
tree.35 <- read.tree(text = TREE)

KrzMCMC = function(Ps, tree, sample.size) {
  # First we calculate the mean a posteriory estimates for the population matrices
  avgPs = llply(alply(apply(Ps, 1:3, mean), 3), function(x) as.matrix(nearPD(x)[[1]]))
  names(avgPs) = dimnames(Ps)[[3]]
  # Then, using AncestralStates, we calculate the within group W-matrix.
  avgP = AncestralStates(tree, avgPs, sample.size)[[tree$edge[1,1]]]
  m = dim(Ps)[3]
  n = dim(Ps)[1]
  # This W-matrix is assumed to be share between all popualtions to form
  # the null hipotesis, and we draw 100 Wishart samples for each population,
  # representing a null distribution of matrices.
  rand.Ps = laply(dimnames(Ps)[[3]],
                  function(x) MonteCarloStat(avgP,
                                             samples[which(dimnames(Ps)[[3]]==x)], 100,
                                             ComparisonFunc=function(x, y) y, cov))
  rand.Ps = aperm(rand.Ps, c(3, 4, 1, 2))
  # We then calculate the shared subspace for the observed and random samples.
  MCMCG.kr.xenartra <- kr.subspace(Ps, vec = rep(n/2, m))
  MCMCG.kr.rand <- kr.subspace(rand.Ps, vec = rep(n/2, m))
  return(list(obs = MCMCG.kr.xenartra,
              rand = MCMCG.kr.rand))
}

load('./maindata.RData')
sample.size = ldply(x, function(x) as.data.frame(x$df))

MCMC.kr = list()

load('./xenartraMCMC.25.samples.Rdata')
MCMC.kr[['25']] = KrzMCMC(Ps, tree.25, sample.size$D25)

load('./xenartraMCMC.28.samples.Rdata')
MCMC.kr[['28']] = KrzMCMC(Ps, tree.25, sample.size$D28)

load('./xenartraMCMC.32.samples.Rdata')
MCMC.kr[['32']] = KrzMCMC(Ps, tree.35, sample.size$D32)

load('./xenartraMCMC.35.samples.Rdata')
MCMC.kr[['35']] = KrzMCMC(Ps, tree.35, sample.size$D35)

dat.krz = ldply(names(MCMC.kr),
                function(x) KrzSubspaceDataFrame (MCMC.kr[[x]][[1]],
                                                  MCMC.kr[[x]][[2]],
                                                  as.numeric(x)))

krz.plot = PlotKrzSubspace(dat.krz)
print(krz.plot)

# divergent.eigen = round(eigen(MCMCG.kr.xenartra$avH)$vectors,3)[,1:10]
# rownames(divergent.eigen) = dimnames(Ps)[[1]]
#
# eigen.correlation = abs(cos(round(apply(MCMCG.kr.xenartra$MCMC.H.theta,1:2,mean),1))[1:10,])
