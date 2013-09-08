# ------------------------------------------------------------------------------
# Function 'seg'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# Last update: 7 September 2013
# Depends: -
# ------------------------------------------------------------------------------
seg <- 
  function(data, nb, tol = .Machine$double.eps) {

  if (ncol(data) > 2)
    warning("'data' has more than two columns; only the first two will be used")
  
  # Duncan and Duncan's index of dissimilarity
  b <- data[,1]/sum(data[,1])     # Blacks
  w <- data[,2]/sum(data[,2])     # Whites
  d <- sum(abs(b-w))/2

  if (!missing(nb)) {     # If information on geographic distance between
                          # spatial units is provided:
    if (!is.matrix(nb))
      stop("'nb' must be a matrix", call. = FALSE)
    else if (nrow(nb) != ncol(nb))
      stop("'nb' must be a square matrix", call. = FALSE)
    else if (nrow(nb) != nrow(data))
      stop("nrow(nb) must match nrow(data)", call. = FALSE)
    
    # Black proportions in census tracts
    z <- data[,1]/(apply(data, 1, sum) + .Machine$double.eps)
    # Additional spatial component value
    spstr <- 0
    nbvec <- as.vector(nb)
    INDEX <- which(nbvec != 0)
    for (i in 1:length(INDEX)) {
      rowID <- INDEX[i] %% nrow(nb)
      colID <- INDEX[i] %/% nrow(nb)
      if (rowID == 0)
        rowID <- nrow(nb)
      else
        colID <- colID + 1
      spstr <- spstr + (abs(z[colID] - z[rowID]) * nbvec[INDEX[i]])
    }
    d <- d - spstr
  }
  
  return(as.vector(d))
}