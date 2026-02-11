setwd("/home/simon/Documents/Projects/dMod/deriv2Tests/wd")
set.seed(1)
m<-4;k<-3;n<-5;l<-2

# [M,K,B] x [K,N]
A <- array(rnorm(m*k*l),c(m,k,l))
B <- matrix(rnorm(k*n),k,n)
da <- dim(A)
db <- dim(B)
C <- dMod:::bmm_lb(A, B, da[3], da[1], da[2], db[2])
