library(Morphometrics)
library(reshape2)
library(gdata)

p = 15
GenerateABCD <- function(p = 15, a_inflation = 10, b_var = 2){
    ret.dim = floor(p/2) -1
    mat_D = RandomMatrix(p)
    eigen_D = eigen(mat_D)
    mat_A = ChangeEigen(mat_D, c(eigen_D$values[1]*a_inflation, eigen_D$values[2:p]))
    mat_B = ChangeEigen(mat_D, sort(abs(rnorm(15, 0, b_var)), decreasing = TRUE))
    eigen_A = eigen(mat_A)
    mat_C = ChangeEigen(mat_A, c(sort(eigen_A$values[1:ret.dim]), eigen_A$values[(ret.dim+1):p]))
    mat_D = RandomMatrix(p, variance = abs(rnorm(15, 0, b_var)))
    matrices = list(A = mat_A, B = mat_B, C = mat_C, D = mat_D)
    RS = melt(RandomSkewers(matrices)[[1]])
    RS = RS[RS$value != 0, ]
    krz = melt(KrzCor(matrices))
    krz = krz[krz$value != 0, ]
    RS$rs = RS$value
    RS$krz = krz$value
    RS$comparison = paste(RS$Var1, RS$Var2, sep = '-')
    RS$value = RS$Var1 = RS$Var2 = NULL
    rownames(RS) = 1:6
    return(RS)
}

ChangeEigen <- function(mat, new_eigen){
    eigen_mat = eigen(mat)
    new_mat = eigen_mat$vectors %*% diag(new_eigen) %*% t(eigen_mat$vectors)
    return(new_mat)
}
