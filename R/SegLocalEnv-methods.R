# ------------------------------------------------------------------------------
# Methods for class 'SegLocalEnv'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# S3 methods
# ------------------------------------------------------------------------------
# "[[.SegLocalEnv" <- function(i, ...) {
#   validObject(x)
#   slotnames <- slotNames(x)
# 
#   if (is.numeric(i)) {
#     i <- as.integer(i)
#     if (i > length(slotnames))
#       chosen <- NULL
#     else {
#       chosen <- slotnames[i]
#       chosen <- paste("x@", chosen, sep = "")
#       chosen <- eval(parse(text = chosen))
#     }
#   }
#   
#   else if (is.character(i)) {
#     chosen <- paste("x@", i, sep = "")
#     chosen <- eval(parse(text = chosen))
#   }
#   
#   else {
#     chosen <- NULL
#   }
# }

print.SegLocalEnv <- function(x, ...) {
  validObject(x)
  cat("Class                 :", class(x), "\n")
  cat("Number of data points :", nrow(x@coords), "\n")
  cat("Number of data columns:", ncol(x@data), "\n")
  cat("Projection            :", x@proj4string@projargs, "\n")
  cat("Slot names            :", slotNames(x), "\n")
}

summary.SegLocalEnv <- function(object, ...) {
  validObject(object)
  msg <- paste("An object of class \"", class(object), "\"\n", sep = "")
  cat(msg)
  cat("Coordinates:\n")
  tmp <- t(apply(object@coords, 2, range))
  rownames(tmp) <- c("x", "y")
  colnames(tmp) <- c("min", "max")
  print(tmp)
  if (is.na(object@proj4string@projargs)) {
    cat("Is projected: FALSE\n")
  } else {
    cat("Is projected: TRUE\n")
    cat("proj4string : [", object@proj4string@projargs, "]\n", sep = "")
  }
    
  cat("\nData values (%):\n")
  tmp <- t(apply(object@data, 1, function(z) z/sum(z))) * 100
  tmp <- apply(tmp, 2, function(z) summary(z, ...))
  colnames(tmp) <- colnames(object@data)
  print(tmp)

  cat("\nLocal environment composition (%):\n")
  tmp <- t(apply(object@env, 1, function(z) z/sum(z))) * 100
  tmp <- apply(tmp, 2, function(z) summary(z, ...))
  colnames(tmp) <- colnames(object@env)
  print(tmp)
}

update.SegLocalEnv <- function(object, coords, data, env, proj4string, ...) {
  validObject(object)
  if (missing(coords))
    coords <- object@coords  
  if (missing(data))
    data <- object@data
  if (missing(env))
    env <- object@env
  if (missing(proj4string))
    proj4string <- object@proj4string
              
  SegLocalEnv(coords, data, env, proj4string)
}

as.list.SegLocalEnv <- function(x, ...) {
  validObject(x)
  list(coords = x@coords, data = x@data, env = x@env, 
       proj4string = x@proj4string)
}

plot.SegLocalEnv <- function(x, which.col, ...) {
  validObject(x)
  xx <- x@coords[,1]
  yy <- x@coords[,2]
  env <- x@env
  
  if (missing(which.col))
    which.col <- 1:ncol(env)
  numPlot <- length(which.col) 
   
  for (i in 1:numPlot) {
    qq <- quantile(env[,i], probs = c(0.2, 0.4, 0.6, 0.8))
    brks <- length(qq) + 1
    size <- rep(brks, nrow(env))
    for (j in 1:brks)
      size[which(env[,i] <= qq[brks - j])] <- brks - j
    plot(x = xx, y = yy, cex = size, ...)
  }
}

points.SegLocalEnv <- function(x, which.col = 1, ...) {
  validObject(x)
  xx <- x@coords[,1]
  yy <- x@coords[,2]
  env <- x@env
  
  i <- which.col[1] 
   
  qq <- quantile(env[,i], probs = c(0.2, 0.4, 0.6, 0.8))
  brks <- length(qq) + 1
  size <- rep(brks, nrow(env))
  for (j in 1:brks)
    size[which(env[,i] <= qq[brks - j])] <- brks - j
  points(x = xx, y = yy, cex = size, ...)
}

# ------------------------------------------------------------------------------
# S4 methods
# ------------------------------------------------------------------------------
setMethod("show", signature(object = "SegLocalEnv"), function(object) {
  validObject(object)
  cat("Class                 :", class(object), "\n")
  cat("Number of data points :", nrow(object@coords), "\n")
  cat("Number of data columns:", ncol(object@data), "\n")
  cat("Projection            :", object@proj4string@projargs, "\n")
  cat("Slot names            :", slotNames(object), "\n")
})

# Coercion methods
#   1) Class 'list' to 'SegLocalEnv', and vice versa
setAs("list", "SegLocalEnv",
      function(from) {
        if (is.null(from$proj4string))
          SegLocalEnv(coords = from$coords, data = from$data, env = from$env)
        else
          SegLocalEnv(coords = from$coords, data = from$data, env = from$env, 
                      proj4string = from$proj4string)
      })
setAs("SegLocalEnv", "list", 
  function(from) {
    validObject(from)
    list(coords = from@coords, data = from@data, env = from@env, 
         proj4string = from@proj4string)
  })
      
#   2) Class 'SpatialPointsDataFrame' to 'SegLocalEnv', and vice versa
setAs("SpatialPointsDataFrame", "SegLocalEnv", 
      function(from) {
        SegLocalEnv(coords = from@coords, data = from@data, env = from@data, 
                    proj4string = from@proj4string)
      })

setAs("SegLocalEnv", "SpatialPointsDataFrame", 
      function(from) {
        validObject(from)
        SpatialPointsDataFrame(coords = from@coords, 
                               data = as.data.frame(from@env),
                               proj4string = from@proj4string)
      })

#   3) Class 'SpatialPolygonsDataFrame' to 'SegLocalEnv'. The other way around
#      is not possible as we have only a set of (x, y) coordinates.
setAs("SpatialPolygonsDataFrame", "SegLocalEnv", 
      function(from) {
        SegLocalEnv(coords = coordinates(from), data = from@data, 
                    env = from@data, proj4string = from@proj4string)
      })
