library(plyr)
library(ggplot2)
library(MCMCglmm)

kr.subspace <- function(Gs, vec = NULL){

  n <- dim(Gs)[[1]]
  m <- dim(Gs)[[3]]
  MCMCsamp <- dim(Gs)[[4]]

  if(is.null(vec)){
    val <- matrix(,n,m)
    for (i in 1:m){
      avG <- apply(Gs, 1:3, mean)
      val[,i] <- round(cumsum(t(eigen(avG[,,i])$values))/sum(eigen(avG[,,i])$values)*100)
    }

    vec <- apply(ifelse(round(val,1) < 90, 1, 0),2,sum)+ 1
  }

  if (dim(Gs)[[1]] != dim(Gs)[[2]]){
    stop("G array must be of order n x n x m x MCMCsamp")
  }

  if(length(vec) != m){stop("vec must have length = m")}
  h <- function (g ,v){
    AA <- array(,c(n, n, m))
    for (k in 1:m){
      g.vec <- eigen(g[,,k])$vectors[,1:(v[k])]
      AA[,,k] <- g.vec %*% t(g.vec)
    }
    H <- apply(AA,1:2,sum)
    list(H = H, AA = AA)
  }
  MCMC.H <- array(,c(n,n,MCMCsamp))
  dimnames(MCMC.H) <- list(dimnames(Gs)[[1]],dimnames(Gs)[[1]],dimnames(Gs)[[4]])
  MCMC.AA <- array(,c(n,n,m,MCMCsamp))
  dimnames(MCMC.AA) <- list(dimnames(Gs)[[1]],dimnames(Gs)[[1]],dimnames(Gs)[[3]],dimnames(Gs)[[4]])
  for (i in 1:MCMCsamp){
    kr <- h(Gs[,,,i], v = vec)
    MCMC.H[,,i] <- kr$H
    MCMC.AA[,,,i] <- kr$AA
  }
  avH <- apply(MCMC.H,1:2,mean)
  rownames(avH) <- dimnames(Gs)[[1]]
  colnames(avH) <- dimnames(Gs)[[1]]
  avAA <- apply(MCMC.AA,1:3,mean)
  dimnames(avAA) <- list(dimnames(Gs)[[1]],dimnames(Gs)[[1]],dimnames(Gs)[[3]])
  avH.vec <- eigen(avH)$vectors
  proj<- function(a,b) t(b) %*% a %*% b
  avH.theta <- matrix(,n,m)
  for (i in 1:n){
    for (i in 1:n){
      avH.theta[i,] <- acos(sqrt(apply(avAA, 3, proj, b = avH.vec[,i]))) * (180/pi)
    }
  }
  MCMC.H.val <- matrix(,MCMCsamp,n)
  colnames(MCMC.H.val) <- paste("h",1:n,sep="")
  for (i in 1:n){
    MCMC.H.val[,i] <- apply(MCMC.H, 3, proj, b = avH.vec[,i])
  }
  MCMC.H.theta <- array(,c(n,m,MCMCsamp))
  rownames(MCMC.H.theta) <- paste("h",1:n,sep="")
  colnames(MCMC.H.theta) <- dimnames(Gs)[[3]]
  for(i in 1:n){
    for(j in 1:MCMCsamp){
      MCMC.H.theta[i,,j] <- acos(sqrt(apply(MCMC.AA[,,,j], 3, proj, b = avH.vec[,i]))) * (180/pi)
    }
  }
  list(avAA = avAA, avH = avH, MCMC.AA = MCMC.AA, MCMC.H = MCMC.H, MCMC.H.val = MCMC.H.val, MCMC.H.theta = MCMC.H.theta)
}

KrzSubspaceDataFrame <- function(MCMCG.kr, MCMCG.kr.rand, n){
  HPD.H.val <- cbind(HPDinterval(as.mcmc(MCMCG.kr$MCMC.H.val)),
                     HPDinterval(as.mcmc(MCMCG.kr.rand$MCMC.H.val)))

  dat.observed = as.data.frame(HPD.H.val[,1:2])
  dat.random = as.data.frame(HPD.H.val[,3:4])
  dat.observed$class = 'observed'
  dat.random$class = 'random'
  dat.observed$PC = dat.random$PC = paste('PC', 1:n, sep = '')
  dat = rbind(dat.observed, dat.random)
  dat = mutate(dat, mean = rowMeans(cbind(lower, upper)))

  orderlist = paste('PC', 1:n, sep = '')
  dat$PC = factor(dat$PC, orderlist)
  dat$n.traits = paste0(n, "D")

  return(dat)
}

PlotKrzSubspace = function(dat){
  plot = ggplot(dat, aes(PC, mean, color = class, group = interaction(PC, class))) + geom_point() +
    geom_errorbar(aes(ymin=lower, ymax=upper)) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab('eigenvectors') + ylab('eigenvalues') + facet_wrap(~ n.traits, ncol = 2, scales = "free_y")
  return(plot)
}
