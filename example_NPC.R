
set.seed(123)
TT <- 10
off_diag <- rnorm(TT, mean=0, sd= 3)
C <- GFT_inverse_mapping(off_diag, tol_value= 10^(-10))
C

library(fBasics)
isPositiveDefinite(C$C)