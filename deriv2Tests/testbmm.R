rm(list = ls(all.names = TRUE))
# Create and set a specific working directory inside your project folder
.workingDir <- file.path(purrr::reduce(1:1, ~dirname(.x), .init = rstudioapi::getSourceEditorContext()$path), "wd")
if (!dir.exists(.workingDir)) dir.create(.workingDir)
setwd(.workingDir)
set.seed(1)
m<-4;k<-3;n<-5;l<-2

# [M,K,B] x [K,N]
A <- array(rnorm(m*k*l),c(m,k,l))
B <- matrix(rnorm(k*n),k,n)
da <- dim(A)
db <- dim(B)
C <- dMod:::bmm_lb(A, B, da[3], da[1], da[2], db[2])

# Comparison:
C_ref <- array(0, c(m, n, l))
for (i in 1:l) C_ref[,,i] <- A[,,i] %*% B

all.equal(C, C_ref)

set.seed(1)
m <- 4; k <- 3; n <- 5; l <- 2

# [M,K,B] x [K,N,B]
A <- array(rnorm(m*k*l), c(m, k, l))
B <- array(rnorm(k*n*l), c(k, n, l))

da <- dim(A)
db <- dim(B)
C <- dMod:::bmm_bb(A, B, da[3], da[1], da[2], db[2])

# Comparison:
C_ref <- array(0, c(m, n, l))
for (i in 1:l) C_ref[,,i] <- A[,,i] %*% B[,,i]

all.equal(C, C_ref)
