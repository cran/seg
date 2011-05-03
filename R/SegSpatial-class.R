# ------------------------------------------------------------------------------
# Class 'SegSpatial'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
setClass("SegSpatial", representation(d = "numeric", r = "numeric", 
                                      h = "numeric", p = "matrix"), 
                       prototype = list(d = numeric(), r = numeric(),
                                        h = numeric(), p = matrix(0, 0, 0)))

setValidity("SegSpatial", function(object) { 
                            if (nrow(object@p) != ncol(object@p))
                              paste("'p' should be a rectangular matrix")
                            else if (!is.numeric(object@p))
                              paste("'p' should be a numeric matrix")
                            else
                              TRUE
                          })
