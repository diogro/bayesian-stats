R.proj <- function(Gs, p = 0.95, vec = 1000){
  if (dim(Gs)[[1]] != dim(Gs)[[2]]){
    stop("G array must be of order n x n x m x MCMCsamp")
  }
  if (is.na(dim(Gs)[4])) {
    stop("There are no MCMCsamples")
  }
  n <- dim(Gs)[[1]]
  m <- dim(Gs)[[3]]
  MCMCsamp <- dim(Gs)[[4]]
  rand.vec <-matrix(,vec,n)
  for (i in 1:vec){
    b <- runif(n,-1,1)
    rand.vec[i,] <- b/(sqrt(sum(b^2)))
  }
  #generate unit length random vectors
  proj<- function(G,b) t(b) %*% G %*% (b)
  #internal fuction to do projection
  G.proj <- array(,c(MCMCsamp,m,vec))
  colnames(G.proj) <- dimnames(Gs)[[3]]
  for (i in 1:vec){
    G.proj[,,i]<- t(apply(Gs,3:4,proj,b = rand.vec[i,]))
  }
  #project each random vector through each MCMC sample of each G
  prs <- cbind(rep(1:m, each = m), 1:m)
  prs.comp <- prs[prs[, 1] < prs[, 2], , drop = FALSE]
  #setting up an index for HPD comparisons
  proj.score <-matrix(,vec,((m^2 - m)/2))
  for (k in 1:vec){
    HPD.int <- HPDinterval(as.mcmc(G.proj[,,k]), prob = p)
    proj.score[k,] <- ifelse(HPD.int[prs.comp[,1],1] > HPD.int[prs.comp[,2],2] | HPD.int[prs.comp[,2],1] > HPD.int[prs.comp[,1],2],1,0)
  }
  #for a given random vector, examine if the HPDintervals of any pair of  G matrices overlap
  vec.score <-cbind(rand.vec,proj.score)
  colnames(vec.score) <- c(1:n, paste(dimnames(Gs)[[3]][prs.comp[, 1]], ".vs.", dimnames(Gs)[[3]][prs.comp[, 2]], sep = ""))
  #collate the random vectors and the outcome of their projection on the G matrices
  sig.vec <- subset(vec.score,rowSums(vec.score[,(n+1):(n+((m^2 - m)/2))]) > 0)
  #collate just the random vectors that resulted signficant differences in variance
  if(dim(sig.vec)[1] <= 1) {warning("There were <= 1 significant vectors, try a larger vec or lower p"); eig.R <- "Na"}
  else{
    eig.R <- eigen(cov(sig.vec[,1:n]))
    rownames(eig.R$vectors) <- dimnames(Gs)[[1]]
    colnames(eig.R$vectors) <- c(paste("e", 1:n, sep = ""))
  }
  #eigen analysis of the R matrix
  list(G.proj = G.proj, vec.score = vec.score, eig.R = eig.R)
}

PlotBayesianRS <- function (MCMC.R.proj, Ps, ncols = 5){
  library(reshape2)

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
  names(dat) = c('Genus', 'interval', 'trait', 'value')
  dat = dcast(dat, Genus+trait~interval)
  names(dat) = c('Genus', 'trait', 'lower', 'upper')
  dat$trait = paste0("PC", dat$trait)
  order.list = paste0("PC", 1:n)
  dat$trait = factor(dat$trait, levels = order.list)
  # levels(dat$pop) <- Gnames
  dat$mean = rowMeans(cbind(dat$upper, dat$lower))
  plot = ggplot(dat, aes(Genus, mean)) + geom_point() +
      geom_errorbar(aes( ymin=lower, ymax=upper)) + theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
     ylab('Genetic Variance') + facet_wrap(~ trait, ncol = ncols, scales = "free_y")
  return(plot)
}
