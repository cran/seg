# ------------------------------------------------------------------------------
# Methods for class 'SegSpatialExt'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# S3 methods
# ------------------------------------------------------------------------------
print.SegSpatialExt <- function(x, digits = getOption("digits"), ...) {
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
    cat("\n")
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

plot.SegSpatialExt <- function(x, smooth = TRUE, which.col, main, ...) {
  validObject(x)
  if (missing(which.col))
    which.col <- 1:ncol(x@localenv@data)
  numPlot <- length(which.col) 
  if (missing(main))
    main <- paste("Data", 1:numPlot)
  else if (length(main) < numPlot)
    main <- rep(main, ceiling(numPlot/length(main)))

  if (!smooth)
    plot(x@localenv, which.col = which.col, main = main, ...)
  else {
    for (i in 1:numPlot)
      plot(x@smooth[[i]], main = main[i], ...)
  }
}

as.list.SegSpatialExt <- function(x, ...) {
  validObject(x)
  list(d = x@d, r = x@r, h = x@h, p = x@p)
}

# ------------------------------------------------------------------------------
# S4 methods
# ------------------------------------------------------------------------------
setMethod("show", signature(object = "SegSpatialExt"), function(object) {
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
    cat("\n")
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
setAs("SegSpatialExt", "list", 
      function(from) {
        validObject(from)
        list(d = from@d, r = from@r, h = from@h, p = from@p)
      })
