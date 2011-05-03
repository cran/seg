# ------------------------------------------------------------------------------
# Function 'getSegLocalEnv'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
getSegLocalEnv <- 
  function(x, data, power = 2, useExp = TRUE, maxdist, sprel, 
           error = .Machine$double.eps) {
  
  # If 'x' is an object of class "Spatial" or one that extends the class
  # (e.g., SpatialPoints), then do the following:
  if (inherits(x, "Spatial")) {
    coords <- try(coordinates(x), silent = TRUE)
    if (class(coords) == "try-error")
      stop("failed to extract coordinates from 'x'", call. = FALSE)
    
    if (missing(data)) {
      # The code below will success if 'x' includes a data frame (e.g., 'x' is
      # a SpatialPointsDataFrame object).
      data <- try(as.matrix(x@data), silent = TRUE) 
      if (class(data) == "try-error")
        stop("'data' is missing, with no default", call. = FALSE)
    }
    removeNA <- which(apply(data, 1, function(z) any(is.na(z))))
    if (length(removeNA) > 0) {
      data <- data[-removeNA,]
      coords <- coords[-removeNA,]
    }

    proj4string <- CRS(as.character(x@proj4string@projargs))
  }

  # If 'x' is an object of class "ppp", then do the following:
  else if (is(x, "ppp")) {
    coords <- cbind(x = x$x, y = x$y)
    proj4string <- CRS(as.character(NA))
    if (missing(data)) {
      if (x$mark != "none")
        stop("'data' is missing, with no default", call. = FALSE)
      else
        data <- as.matrix(x$mark)
    }
  }
  
  # If 'x' is a n * 2 matrix or data frame object, then do the following:
  else if (is.matrix(x) || is.data.frame(x)) {
    coords <- as.matrix(x)
    proj4string <- CRS(as.character(NA))
    if (ncol(coords) != 2 || !is.numeric(coords))
      stop("'x' must be a numeric matrix with two columns", call. = FALSE)
    if (missing(data))
      stop("'data' is missing, with no default", call. = FALSE)
    data <- as.matrix(data)
  }
  
  # If 'x' is not one of the supporting classes:
  else {
    stop("invalid object 'x'", call. = FALSE)
  }

  if (missing(maxdist))
    maxdist <- -1
  else if (!is.numeric(maxdist))
    stop("'maxdist' must be numeric", call. = FALSE)
  else if (maxdist < 0)
    stop("'maxdist' must be greater than or equal to 0", call. = FALSE)
  
  if (missing(sprel)) 
    sprel <- coords
  else {
    if (class(sprel) != "nb" && class(sprel) != "dist")
      stop("invalid object 'sprel'", call. = FALSE)
  }
  
  env <- .SEGLOCALENV(sprel, data, power, useExp, maxdist, error)
  SegLocalEnv(coords, data, env, proj4string)
}

# ------------------------------------------------------------------------------
# Internal function '.SEGLOCALENV'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
.SEGLOCALENV <- function(x, data, power, useExp, maxdist, error) {

  if (ncol(data) < 2 || !is.numeric(data))
    stop("'data' must be a numeric matrix with at least two columns", 
         call. = FALSE)

  if (!is.numeric(power))
    stop("'power' must be numeric", call. = FALSE)
  if (!is.logical(useExp))
    stop("'useExp' must be logical", call. = FALSE)
  if (!is.numeric(error))
    stop("'error' must be numeric", call. = FALSE)
  
  if (class(x) == "nb") {
    # tt <- require(spdep)
    # if (!tt)
    #   stop("package 'spdep' cannot be loaded", call. = FALSE)
    xmat <- nb2mat(x, style = "W")
    if (nrow(xmat) != nrow(data)) 
      stop("'data' must have the same number of rows as 'x'", call. = FALSE)
    env <- xmat %*% data
  } else if (class(x) == "dist") {
    x <- as.matrix(x)
    if (nrow(x) != nrow(data)) 
      stop("'data' must have the same number of rows as 'x'", call. = FALSE)
    env <- matrix(nrow = nrow(data), ncol = ncol(data))
    for (i in 1:nrow(data)) {
      if (useExp)
        weight <- exp(power * x[i,] * -1)
      else
        weight <- 1/(x[i,] + error)^power
      if (maxdist >= 0)
        weight[which(x[i,] > maxdist)] <- 0
      env[i,] <- apply(data, 2, function(z) sum(z * weight)/sum(weight))
    }
  } else {
    if (nrow(x) != nrow(data))
      stop("'data' must have the same number of rows as 'x'", call. = FALSE)
    xval <- x[,1]; yval <- x[,2]
    dim <- ncol(data); data <- as.vector(data)
    env <- .Call("envconstruct", xval, yval, data, as.integer(dim), power, 
                 as.integer(useExp), maxdist, error)
    # --------------------------------------------------------------------------
    # R version 'envconstruct()'
    # --------------------------------------------------------------------------
	  # envconstruct <- function(x, y, v, d) {
    #   n <- length(v)
    #   env <- rep(0, n)
    #   
    #   for (i in 1:n) {
    #     m <- 0
    #     for (j in 1:n) {
    #       dx <- x[i] - x[j]
    #       dy <- y[i] - y[j]
    #       if (dx <= d && dy <= d) {
    #         if (dx^2 + dy^2 <= d^2) {
    #           env[i] <- env[i] + v[j]
    #           m <- m + 1
    #         }
    #       }
    #     }
    #     if (m > 1)
    #       env[i] <- env[i] / m
    #   }
    #   return(env)
    # }
    # --------------------------------------------------------------------------
  }
  invisible(env)
}
