# Generate eqnlist from the constructor
S <- matrix(c(-1, 1, 1, -1),
            nrow = 2, ncol = 2,
            dimnames = list(NULL, c("A", "B")))

rates <- c("k1*A", "k2*B")
description <- c("forward", "backward")

f <- eqnlist(smatrix = S, rates = rates, description = description)
print(f)

# Convert to data.frame
fdata <- as.data.frame(f)
print(fdata)

# Legacy path: flat `volumes` gets auto-translated into compartments c1, c2
f <- as.eqnlist(fdata, volumes = c(A = "Vcyt", B = "Vnuc"))
print(f)
print(as.eqnvec(f))
print(as.eqnvec(f, type = "amount"))

# First-class compartments: name each compartment and assign states explicitly.
f2 <- eqnlist(
  smatrix = S, rates = rates, description = description,
  compartments = list(cyt = "Vcyt", nuc = "Vnuc"),
  compartmentOf = c(A = "cyt", B = "nuc")
)
print(f2)
print(as.eqnvec(f2))
