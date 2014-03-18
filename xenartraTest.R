load("./maindata.RData")

getModel <- function(otu){
  models = otu$model[!is.na(otu$model)]
  models = models[[1]]
  return(as.character(models$call$formula)[[3]])
}

models = llply(x, getModel)

for(i in 1:length(x)){
  x[[i]]$data <- cbind(x[[i]]$raw, x[[i]]$categorical)
  x[[i]]$fixed <- models[[i]]
}