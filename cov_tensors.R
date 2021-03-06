library(plyr)
library(MCMCglmm)
library(gdata)
library(matrixcalc)
library(ggplot2)
library(Morphometrics)
library(reshape2)

BuildSigmaMat <- function(avG, n){
  neigten <- n*(n+1)/2
  m = dim(avG)[[3]]
  Smat <- matrix(nrow=neigten, ncol=neigten)
  varmat <- matrix(nrow=m, ncol=n)
  covmat <- matrix(nrow=m, ncol=(n*n-n)/2)
  for(i in 1:m){
    varmat[i,] = diag(avG[,,i])
    covmat[i,] = lowerTriangle(avG[,,i])
  }
  Smat[1:n,1:n] <- cov(varmat,varmat)
  #Fill the upper left quadrant of S
  Smat[(n+1):neigten,(n+1):neigten] <-2*cov(covmat,covmat)
  #Fill the lower right quadrant of S
  Smat[1:n,(n+1):neigten] <- sqrt(2) * cov(varmat,covmat)
  #Fill the upper right quadrant of S
  Smat[(n+1):neigten,1:n] <- sqrt(2) * cov(covmat,varmat)
  return(Smat)
}

covtensor <- function(Gs){
  dimsG <- dim(Gs)
  n = dimsG[1]
  m = dimsG[3]
  MCMCsamp = dimsG[4]
  avG <- apply(Gs, 1:3, mean)
  #Calculate avG
  neigten <- n*(n+1)/2
  #Number of eigentensors
  Smat = BuildSigmaMat(avG, n)
  #Fill the lower left quadrant of S
  Seigvals <- eigen(Smat)$values
  #Eigenvalues of S
  Seigvecs <- eigen(Smat)$vectors
  #Eigenvectors of S
  eTmat <- array(,c(n, n, neigten))
  for (i in 1:neigten){
    emat <- matrix(0,nrow=n,ncol=n)
    lowerTriangle(emat) <- 1/sqrt(2)*Seigvecs[(n+1):neigten,i]
    emat <- emat + t(emat)
    diag(emat) <- Seigvecs[1:n,i]
    eTmat[,,i] <- emat
    #Using the eigenvectors of S construct the second-order eigentensors
  }
  eTeigen <- array(,c(n+1, n, neigten))
  for (i in 1:neigten){
    eTeigen[1,,i] <- t(eigen(eTmat[,,i])$values)
    #Eigenvalues of the ith eigentensor
    eTeigen[2:(n+1),,i] <- eigen(eTmat[,,i])$vectors
    #Eigenvectors of the ith eigentensor
  }
  avG.coord <- array(,c(m, neigten, 1))
  for (i in 1:neigten){
    avG.coord[,i,] <- apply(avG, 3, frobenius.prod, y = eTmat[,,i])
    #Coordinates of the jth avG for the ith eigentensor
  }
  if(is.na(dim(Gs)[4])) {warning("There are no MCMCsamples"); MCMCSeigvals <- "NA";  MCMCG.coord <- "NA"}  else {
    MCMCSeigvals <- matrix(,MCMCsamp,neigten)
    for (k in 1:MCMCsamp)
      for (l in 1:neigten)
        MCMCSeigvals[k,l] <- t(Seigvecs[,l]) %*% BuildSigmaMat(Gs[,,,k], n) %*% Seigvecs[,l]

    MCMCG.coord <- array(,c(m, neigten, MCMCsamp))
    for (i in 1:neigten){
      MCMCG.coord[,i,] <- apply(Gs, 3:4, frobenius.prod, y = eTmat[,,i])
      #Coordinates of the kth MCMC sample of the jth G for the ith eigentensor
    }
  }
  tensor.summary <- data.frame(rep(Seigvals,each=n),rep(propeigval<-Seigvals/sum(Seigvals),each=n),t(data.frame(eTeigen)))
  if(!exists("traitnames")) traitnames <- letters[1:n]
  colnames(tensor.summary) <- c("eigenvalue1","proportion","eigenvalue2",traitnames)
  ordered.tensor.summary <- tensor.summary[order(tensor.summary$eigenvalue1,abs(tensor.summary$eigenvalue2),decreasing=T),]
  rownames(ordered.tensor.summary) <- paste(paste("tensor",rep(1:neigten, each=n),sep=""),paste("e",rep(1:n,neigten),sep=""))
  list(Smat = Smat, eTmat = eTmat, ordered.tensor.summary = ordered.tensor.summary, avG.coord = avG.coord, MCMCSeigvals = MCMCSeigvals, MCMCG.coord = MCMCG.coord)
}
#END

CovTensorDataFrame <- function(MCMC.covtensor, MCMC.covtensor.rand, n){

  m = sum(!is.na(sample.size[paste('D',n, sep='')]))
  nnonzero <- min(n*(n+1)/2,m-1)

  HPD.eT.val <- cbind(HPDinterval(as.mcmc(MCMC.covtensor$MCMCSeigvals[,1:nnonzero]), prob=0.95),
                      HPDinterval(as.mcmc(MCMC.covtensor.rand$MCMCSeigvals[,1:nnonzero]), prob=0.95))

  dat.observed = as.data.frame(HPD.eT.val[,1:2])
  dat.random = as.data.frame(HPD.eT.val[,3:4])
  dat.observed$class = 'observed'
  dat.random$class = 'random'
  dat.observed$PC = dat.random$PC = paste('PC', 1:nnonzero, sep = '')
  dat = rbind(dat.observed, dat.random)
  dat = mutate(dat, mean = rowMeans(cbind(lower, upper)))

  orderlist = paste('PC', 1:nnonzero, sep = '')
  dat$PC = factor(dat$PC, orderlist)
  dat$n.traits = paste0(n, "D")
  return(dat)
}

PlotCovTensor = function(dat){
  plot = ggplot(dat, aes(PC, mean, color = class, group = PC)) + scale_y_log10() +
    geom_point() + geom_errorbar(aes(ymin=lower, ymax=upper)) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab('eigenvectors') + ylab('eigenvalues') + facet_wrap(~ n.traits, ncol = 2, scales = "free_y")
  return(plot)
}
