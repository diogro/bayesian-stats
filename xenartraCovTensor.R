source('./cov_tensors.R')
source('xenartraSamples.R')


library(Morphometrics)

load('./maindata.RData')
sample.size = ldply(x, function(x) as.data.frame(x$df))

covTensor = function(Ps, tree, sample.size) {
  m = dim(Ps)[3]
  n = dim(Ps)[1]
  SamplePop = function(pop_Ps, n.ind, otu){
    pop = adply(pop_Ps, 3, function(x) mvtnorm::rmvnorm(n = n.ind,
                                                        mean = rep(0, dim(x)[2]),
                                                        sigma = x))
    pop$'.id' = otu
    return(pop)
  }
  pop_Ps = ldply(dimnames(Ps)[[3]], function(otu) SamplePop(Ps[,,which(dimnames(Ps)[[3]]==otu),],
                                                            sample.size[which(dimnames(Ps)[[3]]==otu)],
                                                            otu))
  shuffle = sample(dim(pop_Ps)[1])
  pop_Ps$'.id' = pop_Ps$'.id'[shuffle]
  pop_Ps$'X1' = pop_Ps$'X1'[shuffle]
  rand.Ps = daply(pop_Ps, .(X1, .id), function(x) cov(x[-c(1, dim(x)[2])]))
  rand.Ps = aperm(rand.Ps, c(3, 4, 2, 1))

  MCMC.covtensor <- covtensor(Ps)
  MCMC.covtensor.rand <- covtensor(rand.Ps)
  return(list(obs = MCMC.covtensor,
              rand = MCMC.covtensor.rand))
}

#MCMC.covTensor = list()
#MCMC.covTensor[['25']] = covTensor(Ps[['25']], tree.25, sample.size$D25)
#MCMC.covTensor[['28']] = covTensor(Ps[['28']], tree.25, sample.size$D28)
#MCMC.covTensor[['32']] = covTensor(Ps[['32']], tree.35, sample.size$D32)
#MCMC.covTensor[['35']] = covTensor(Ps[['35']], tree.35, sample.size$D35)
#save(MCMC.covTensor, file = 'xenartraCovTensor.Rdata')
load('xenartraCovTensor.Rdata')

dat.cov.tensor = ldply(names(MCMC.covTensor),
                       function(x) CovTensorDataFrame (MCMC.covTensor[[x]][[1]],
                                                       MCMC.covTensor[[x]][[2]],
                                                       as.numeric(x)))

cov.tensor.plot = PlotCovTensor(dat.cov.tensor)
ggsave("~/Desktop/cov_tensor_plot.tiff")

m = dim(Ps)[3]
n = dim(Ps)[1]
nnonzero <- min(n*(n+1)/2,m-1)

HPD.eT.val <- cbind(HPDinterval(as.mcmc(MCMC.covTensor[[3]][[1]]$MCMCSeigvals[,1:nnonzero]), prob=0.95),
                    HPDinterval(as.mcmc(MCMC.covTensor[[3]][[2]]$MCMCSeigvals[,1:nnonzero]), prob=0.95))
round(HPD.eT.val,3)
round(MCMC.covTensor[[3]][[1]]$ordered.tensor.summary[1:n*2,2:9],3)
