# ------------------------------------------------------------------------------
# Class 'SegLocalEnv'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
setClass("SegLocalEnv", representation(coords = "matrix", data = "matrix", 
                                       env = "matrix", proj4string = "CRS"))

setValidity("SegLocalEnv", 
  function(object) {
    if (ncol(object@coords) != 2 || !is.numeric(object@coords))
      paste("'coords' must be a numeric matrix with two columns, x and y")
    else if (ncol(object@data) < 2 || !is.numeric(object@data))
      paste("'data' must be a numeric matrix with at least two columns")
    else if (nrow(object@data) != nrow(object@coords))
      paste("'data' must have the same number of rows as 'coords'")
    else if (any(dim(object@env) != dim(object@data)))
      paste("'env' must be a matrix with the same dimensions as 'data'")
    else if (!is.numeric(object@env))
      paste("'env' must be a numeric matrix")
    else if (class(object@proj4string) != "CRS")
      paste("'proj4string' is not a valid CRS object")
    else
      TRUE
  })
