source('./xenartraSamples.R')
source('./KrzSubspace.R')

avgP = apply(Ps, 1:2, mean)
m = dim(Ps)[3]
n = dim(Ps)[1]
rand.Ps = laply(dimnames(Ps)[[3]], 
                function(x) MonteCarloStat(avgP, n+1, 100, 
                                           ComparisonFunc=function(x, y) y, cov))
rand.Ps = aperm(rand.Ps, c(3, 4, 1, 2))

MCMCG.kr.xenartra <- kr.subspace(Ps, vec = rep(n/2, m))
MCMCG.kr.rand <- kr.subspace(rand.Ps, vec = rep(n/2, m))

PlotKrzSubspace <- function(MCMCG.kr, MCMCG.kr.rand){
  HPD.H.val <- cbind(HPDinterval(as.mcmc(MCMCG.kr$MCMC.H.val)), HPDinterval(as.mcmc(MCMCG.kr.rand$MCMC.H.val)))
  
  dat.observed = as.data.frame(HPD.H.val[,1:2])
  dat.random = as.data.frame(HPD.H.val[,3:4])
  dat.observed$class = 'observed'
  dat.random$class = 'random'
  dat.observed$PC = dat.random$PC = paste('PC', 1:num.traits, sep = '')
  dat = rbind(dat.observed, dat.random)
  dat = mutate(dat, mean = rowMeans(cbind(lower, upper)))
  
  orderlist = paste('PC', 1:num.traits, sep = '')
  dat$PC = factor(dat$PC, orderlist)
  plot = ggplot(dat, aes(x = reorder(PC, orderlist), mean, color = class, group = PC)) + geom_point() + 
    geom_errorbar(aes(ymin=lower, ymax=upper)) + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    xlab('eigenvectors') + ylab('eigenvalues') 
  return(plot)
}

plot = PlotKrzSubspace(MCMCG.kr.xenartra, MCMCG.kr.rand)
divergent.eigen = round(eigen(MCMCG.kr.xenartra$avH)$vectors,3)[,1:5]
rownames(divergent.eigen) = dimnames(Ps)[[1]]

eigen.correlation = abs(cos(round(apply(MCMCG.kr$MCMC.H.theta,1:2,mean),1))[1:5,])