source('./xenartraSamples.R')
source('./b_random_skewers.R')

PlotBayesianRS <- function (MCMC.R.proj, Ps){
  proj<- function(G,b) t(b) %*% G %*% (b)
  #Function to do projection
  n <- dim(Ps)[[1]]
  m <- dim(Ps)[[3]]
  MCMCsamp <- dim(Ps)[[4]]
  R.vec.proj <- array(,c(MCMCsamp,m,n))
  for (i in 1:n){
    R.vec.proj[,,i] <- t(apply(Ps, 3:4, proj, b = MCMC.R.proj$eig.R$vectors[,i]))
  }
  #Genetic variance in each population in the direction of the eigenvectors of R
  
  HPD.R.vec.proj <- array(,c(m,2,n))
  for (i in 1:n){
    HPD.R.vec.proj[,,i] <- HPDinterval(as.mcmc(R.vec.proj[,,i]), prob = 0.95)    
  }
  dimnames(HPD.R.vec.proj)[[1]] = dimnames(Ps)[[3]]
  #HPD intervals for the genetic variance in each population in the direction of the eigenvectors of R
  
  dat = adply(HPD.R.vec.proj, 1:3)
  names(dat) = c('pop', 'interval', 'trait', 'value')
  dat = dcast(dat, pop+trait~interval)
  names(dat) = c('pop', 'trait', 'lower', 'upper')
  # levels(dat$pop) <- Gnames
  dat$mean = rowMeans(cbind(dat$upper, dat$lower))
  plot.func = function(eigen.number) {
    label = paste('Pc', eigen.number)
    plot = ggplot(subset(dat, trait==eigen.number), aes(pop, mean)) + geom_point() + 
      geom_errorbar(aes( ymin=lower, ymax=upper)) + theme_bw() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      xlab(label) + ylab('Genetic Variance') 
    return(plot)
  }
  plot.list = alply(1:n, 1, function(x) plot.func(x))
  return(plot.list)
}

MCMC.R.proj = R.proj(Ps)
plots = PlotBayesianRS(MCMC.R.proj, Ps)