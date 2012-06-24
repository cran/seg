# ------------------------------------------------------------------------------
# seg()
# ------------------------------------------------------------------------------
seg <- 
  function(data, nb, verbose = FALSE) {

  if (ncol(data) > 2)
    warning("'x' has more than two columns; only the first two will be used")
  b <- data[,1]/sum(data[,1])
  w <- data[,2]/sum(data[,2])
  z <- data[,1]/(apply(data, 1, sum) + .Machine$double.eps)
  d <- sum(abs(b-w))/2
  
  if (!missing(nb)) {  
    if (isSymmetric(nb)) {
      subsets <- t(combn(1:nrow(nb), 2))
      dd <- as.numeric(as.dist(nb))
    } else {
      subsets <- expand.grid(1:nrow(nb), 1:nrow(nb))
      dd <- as.numeric(nb)
    }
    INDEX <- which(dd != 0)
    subsets <- subsets[INDEX,]
    dd <- dd[INDEX]
    # The code below should work fine but slow!
    #
    # dd <- numeric()
    # for (i in 1:nrow(subsets)) {
    #   rr <- subsets[i, 1]
    #   cc <- subsets[i, 2]
    #   dd <- append(dd, nb[rr, cc])
    # }

    dm <- sum(abs(z[subsets[,1]] - z[subsets[,2]]) * dd)
    d <- d - dm
  }
  
  return(d)
}

