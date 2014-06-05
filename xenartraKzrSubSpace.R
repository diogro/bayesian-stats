source('./KrzSubspace.R')
library(Morphometrics)

KrzMCMC = function(Ps, sample.size) {
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

load('./Rdatas/maindata.RData')
sample.size = ldply(x, function(x) as.data.frame(x$df))

MCMC.kr = c('25', '28', '32', '35')
MCMC.kr = alply(MCMC.kr, 1, function(traits) KrzMCMC(Ps[[traits]], sample.size[[paste0('D', traits)]]), .parallel = TRUE)
names(MCMC.kr) <- c('25', '28', '32', '35')

dat.krz = ldply(names(MCMC.kr),
                function(x) KrzSubspaceDataFrame (MCMC.kr[[x]][[1]],
                                                  MCMC.kr[[x]][[2]],
                                                  as.numeric(x)))

krz.plot = PlotKrzSubspace(dat.krz)
print(krz.plot)
ggsave("~/Desktop/krz_projection.tiff", krz.plot, height = 20, width = 30, units = "cm")

divEigen <- function(x, vec=1){
 divergent.eigen = eigen(x$obs$avH)$vectors[,1:vec]
 rownames(divergent.eigen) = dimnames(Ps)[[1]]
 abs(cos(round(apply(x$obs$MCMC.H.theta,1:2,mean),1))[1:vec,])
}
lapply(MCMC.kr, divEigen, 5)
