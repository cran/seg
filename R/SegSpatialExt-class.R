# ------------------------------------------------------------------------------
# Class 'SegSpatialExt'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
setClass("SegSpatialExt", contains = "SegSpatial", 
                          representation(nrow = "numeric", ncol = "numeric", 
                                         method = "character", 
                                         sigma = "numeric",
                                         smooth = "list", 
                                         localenv = "SegLocalEnv"), 
                          prototype = list(nrow = numeric(), ncol = numeric(),
                                           method = character(),
                                           sigma = numeric(), smooth = list()))

setValidity("SegSpatialExt", 
            function(object) { 
              if (!validObject(object@localenv))
                paste("'localenv' is not a valid object")
              else
                TRUE
            })
