source('./KrzSubspace.R')
library(Morphometrics)

KrzMCMC = function(Ps, tree, sample.size) {
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
    # We then calculate the shared subspace for the observed and random samples.
    MCMCG.kr.xenartra <- kr.subspace(Ps, vec = rep(n/2, m))
    MCMCG.kr.rand <- kr.subspace(rand.Ps, vec = rep(n/2, m))
    return(list(obs = MCMCG.kr.xenartra,
                rand = MCMCG.kr.rand))
}

load('./maindata.RData')
sample.size = ldply(x, function(x) as.data.frame(x$df))

MCMC.kr = list()

MCMC.kr[['25']] = KrzMCMC(Ps[['25']], tree.25, sample.size$D25)
MCMC.kr[['28']] = KrzMCMC(Ps[['28']], tree.25, sample.size$D28)
MCMC.kr[['32']] = KrzMCMC(Ps[['32']], tree.35, sample.size$D32)
MCMC.kr[['35']] = KrzMCMC(Ps[['35']], tree.35, sample.size$D35)

dat.krz = ldply(names(MCMC.kr),
                function(x) KrzSubspaceDataFrame (MCMC.kr[[x]][[1]],
                                                  MCMC.kr[[x]][[2]],
                                                  as.numeric(x)))

krz.plot = PlotKrzSubspace(dat.krz)
print(krz.plot)
ggsave("~/Desktop/krz_projection.png", height = 20, width = 30, units = "cm")

# divergent.eigen = round(eigen(MCMCG.kr.xenartra$avH)$vectors,3)[,1:10]
# rownames(divergent.eigen) = dimnames(Ps)[[1]]
#
# eigen.correlation = abs(cos(round(apply(MCMCG.kr.xenartra$MCMC.H.theta,1:2,mean),1))[1:10,])
