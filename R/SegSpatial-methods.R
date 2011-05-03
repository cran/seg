# ------------------------------------------------------------------------------
# Methods for class 'SegSpatial'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# S3 methods
# ------------------------------------------------------------------------------
print.SegSpatial <- function(x, digits = getOption("digits"), ...) {
  validObject(x)
  cat("\n\tReardon and O'Sullivan's spatial segregation measures\n\n")

  cat("Dissimilarity (D)     : ")
  if (length(x@d) > 0)
    cat(round(x@d, digits), "\n")
  else
    cat("-\n")

  cat("Relative diversity (R): ")
  if (length(x@r) > 0)
    cat(round(x@r, digits), "\n")
  else
    cat("-\n")

  cat("Information theory (H): ")            
  if (length(x@h) > 0)
    cat(round(x@h, digits), "\n")
  else
    cat("-\n")

  cat("Exposure/Isolation (P): ")
  if (length(x@p) > 0) {
    cat("\n\n")
    if (is.null(colnames(x@p)))
      colnames(x@p) <- paste("Group", 1:ncol(x@p))
    if (is.null(rownames(x@p)))
      rownames(x@p) <- paste("Group", 1:nrow(x@p))
    print(x@p, digits, ...)
    cat("--\n")
    cat("The exposure/isolation matrix should be read horizontally.\n")
    cat("Read 'help(spseg)' for more details.\n")
  } else {
    cat("-\n")
  }
}

as.list.SegSpatial <- function(x, ...) {
  validObject(x)
  list(d = x@d, r = x@r, h = x@h, p = x@p)
}

# ------------------------------------------------------------------------------
# S4 methods
# ------------------------------------------------------------------------------
setMethod("show", signature(object = "SegSpatial"), function(object) {
  validObject(object)
  cat("\n\tReardon and O'Sullivan's spatial segregation measures\n\n")

  cat("Dissimilarity (D)     : ")
  if (length(object@d) > 0)
    cat(round(object@d, 4), "\n")
  else
    cat("-\n")

  cat("Relative diversity (R): ")
  if (length(object@r) > 0)
    cat(round(object@r, 4), "\n")
  else
    cat("-\n")

  cat("Information theory (H): ")            
  if (length(object@h) > 0)
    cat(round(object@h, 4), "\n")
  else
    cat("-\n")

  cat("Exposure/Isolation (P): ")
  if (length(object@p) > 0) {
    cat("\n\n")
    if (is.null(colnames(object@p)))
      colnames(object@p) <- paste("Group", 1:ncol(object@p))
    if (is.null(rownames(object@p)))
      rownames(object@p) <- paste("Group", 1:nrow(object@p))
    print(object@p)
    cat("--\n")
    cat("The exposure/isolation matrix should be read horizontally.\n")
    cat("Read 'help(spseg)' for more details.\n")
  } else {
    cat("-\n")
  }
})

# Coercion methods
setAs("list", "SegSpatial",
      function(from) {
        tmp <- new("SegSpatial")
        if (!is.null(from$d)) {
          if (!is.numeric(from$d))
            stop("invalid object at 'd' - must be numeric", call. = FALSE)
          if (length(from$d) == 0)
            tmp@d <- from$d
          else
            tmp@d <- from$d[1]
        }
        
        if (!is.null(from$r)) {
          if (!is.numeric(from$r))
            stop("invalid object at 'r' - must be numeric", call. = FALSE)
          if (length(from$r) == 0)
            tmp@r <- from$r
          else
            tmp@r <- from$r[1]
        }
        
        if (!is.null(from$h)) {
          if (!is.numeric(from$h))
            stop("invalid object at 'h' - must be numeric", call. = FALSE)
          if (length(from$h) == 0)
            tmp@h <- from$h
          else
            tmp@h <- from$h[1]
        }
        
        if (!is.null(from$p)) {
          if (!is.matrix(from$p))
            stop("invalid object at 'p' - must be a matrix", call. = FALSE)
          else if (!is.numeric(from$p))
            stop("invalid object at 'p' - must be a numeric matrix", 
                 call. = FALSE)
          else if (nrow(from$p) != ncol(from$p))
            stop("invalid object at 'p' - must be a rectangular matrix", 
                 call. = FALSE)
          else
            tmp@p <- from$p
        }
        
        return(tmp)
      })

setAs("SegSpatial", "list", 
      function(from) {
        validObject(x)
        list(d = x@d, r = x@r, h = x@h, p = x@p)
      })
