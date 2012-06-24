# ------------------------------------------------------------------------------
# whiteseg()
# ------------------------------------------------------------------------------
whiteseg <- 
  function(x, data, nb, method = "euclidean", p = 2, fun, verbose = FALSE) {

  tmp <- .SEGDATA(x, data, verbose)
  coords <- tmp$coords; data <- tmp$data;
  rm(tmp)

  # Verify 'coords' and 'data'
  if (ncol(data) < 2 || !is.numeric(data))
    stop("'data' must be a numeric matrix with at least two columns", 
         call. = FALSE)
  else if (nrow(data) != nrow(coords))
    stop("'data' must have the same number of rows as 'x'", call. = FALSE)

  # If the user did not specify 'dist', do the followings:
  if (missing(nb))
    nb <- as.matrix(dist(x = coords, method = method, p = p))
  
  # If 'fun' is not given, use the default option (i.e., negative exponential):
  if (missing(fun))
    fun <- function(z) { exp(-z) }
  
  # If the distance matrix is symmetric and if all the left-to-right diagonals
  # are 0, we can save some time
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
  
  dfun <- fun(dd)
  nTracts <- apply(data, 1, sum)
  nGroups <- apply(data, 2, sum)

  pA <- list()
  pB <- sum(nTracts[subsets[,1]] * nTracts[subsets[,2]] * dfun) / sum(nGroups)
  for (i in 1:ncol(data))
    pA[[i]] <- sum(data[subsets[,1], i] * data[subsets[,2], i] * dfun) / 
               nGroups[i]
  sum(unlist(pA)) / pB
}
