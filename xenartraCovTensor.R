source('./xenartraSamples.R')
source('./cov_tensors.R')

avgP = apply(Ps, 1:2, mean)
m = dim(Ps)[3]
n = dim(Ps)[1]
rand.Ps = laply(dimnames(Ps)[[3]], 
                function(x) MonteCarloStat(avgP, n+1, 100,
                                           ComparisonFunc=function(x, y) y, cov))
rand.Ps = aperm(rand.Ps, c(3, 4, 1, 2))

MCMC.covtensor <- covtensor(Ps)
MCMC.covtensor.rand <- covtensor(rand.Ps)
save(MCMC.covtensor, MCMC.covtensor.rand, file = 'covtensor-35.RData')

nnonzero <- min(n*(n+1)/2,m-1)

HPD.eT.val <- cbind(HPDinterval(as.mcmc(MCMC.covtensor$MCMCSeigvals[,1:nnonzero]), prob=0.95), 
                    HPDinterval(as.mcmc(MCMC.covtensor.rand$MCMCSeigvals[,1:nnonzero]), prob=0.95))
round(HPD.eT.val,3)
round(MCMC.covtensor$ordered.tensor.summary[1:16,2:9],3)par(mfrow=c(1,1))

PlotCovTensor <- function(MCMC.covtensor, MCMC.covtensor.rand){
  
  HPD.eT.val <- cbind(HPDinterval(as.mcmc(MCMC.covtensor$MCMCSeigvals[,1:nnonzero]), prob=0.95), 
                      HPDinterval(as.mcmc(MCMC.covtensor.rand$MCMCSeigvals[,1:nnonzero]), prob=0.95))
  
  nnonzero <- min(n*(n+1)/2,m-1)
  
  dat.observed = as.data.frame(HPD.eT.val[,1:2])
  dat.random = as.data.frame(HPD.eT.val[,3:4])
  dat.observed$class = 'observed'
  dat.random$class = 'random'
  dat.observed$PC = dat.random$PC = paste('PC', 1:nnonzero, sep = '')
  dat = rbind(dat.observed, dat.random)
  dat = mutate(dat, mean = rowMeans(cbind(lower, upper)))
  
  orderlist = paste('PC', 1:num.traits, sep = '')
  dat$PC = factor(dat$PC, orderlist)
  plot = ggplot(dat, aes(x = reorder(PC, orderlist), mean, color = class, group = PC)) + geom_point() + 
    geom_errorbar(aes(ymin=lower, ymax=upper)) + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    xlab('eigenvectors') + ylab('eigenvalues') 
  plot2 = ggplot(dplyr::filter(dat, PC != 'PC1'), aes(x = reorder(PC, orderlist), mean, color = class, group = PC)) + geom_point() + 
    geom_errorbar(aes(ymin=lower, ymax=upper)) + theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    xlab('eigenvectors') + ylab('eigenvalues') 
  return(list(all.PC = plot, exclude.PC1 = plot2))
}

PlotCovTensor(MCMC.covtensor, MCMC.covtensor.rand)
