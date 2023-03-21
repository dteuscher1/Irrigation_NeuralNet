##
## Create basis functions for spatial NN
##
library(tidyverse)

locs <- read.csv("./x_matrix.csv")
D <- fields::rdist(locs)
A <- 1/D
diag(A) <- 0
A[1:10, 1:10]
Po <- diag(nrow(locs)) - matrix(1, nrow=nrow(A), ncol=ncol(A))/nrow(A)
M <- Po%*%A%*%Po %>% eigen()
pos <- which(M$values>0)
kp.prop <- .5
kp.num <- which((cumsum(M$values[pos])/sum(M$values[pos]))<=kp.prop)
B <- M$vectors[,kp.num]

## Plot
fields::quilt.plot(locs[,1], locs[,2], B[,12])

## Write it out
write.csv(x=B, file="./SpatialBasisFunctions.csv", 
          quote=FALSE, row.names=FALSE)
