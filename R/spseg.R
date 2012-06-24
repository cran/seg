# ------------------------------------------------------------------------------
# Function 'spseg'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
spseg <- 
  function(x, data, method = "all", smoothing = "none", nrow = 100, ncol = 100, 
           window, sigma, useC = TRUE, negative.rm = FALSE, 
           tol = .Machine$double.eps, verbose = FALSE, ...) {
  excTime <- Sys.time()  
  tmp <- .SEGDATA(x, data, verbose)
  coords <- tmp$coords; data <- tmp$data;
  rm(tmp)
  smoothed.df <- list()

  # Verify 'coords' and 'data' and match 'method' and 'smoothing' arguments 
  # against candidate values.
  if (ncol(data) < 2 || !is.numeric(data))
    stop("'data' must be a numeric matrix with at least two columns", 
         call. = FALSE)
  else if (nrow(data) != nrow(coords))
    stop("'data' must have the same number of rows as 'x'", call. = FALSE)

  method <- match.arg(method, c("exposure", "information", "diversity", 
                                "dissimilarity", "all"), several.ok = TRUE)
  smoothing <- 
    match.arg(smoothing, c("none", "kernel", "equal"), several.ok = FALSE)

  ##############################################################################
  #
  # Step 1. Estimate the data surface
  #
  ##############################################################################  
  if (verbose) {
    begTime <- Sys.time()
    cat("Create a data surface of 'x' ...\n")
  }

  # ----------------------------------------------------------------------------
  # (1) Equal smoothing
  # ----------------------------------------------------------------------------
  if (smoothing == "equal") {
    if (!inherits(x, "SpatialPolygons")) {
      msg <- paste("'x' should be a SpatialPolygons object to", 
                   "use the \"equal\" smoothing option")
      stop(msg, call. = FALSE)
    }

    xmn <- min(coords[,1]); xmx <- max(coords[,1])
    ymn <- min(coords[,2]); ymx <- max(coords[,2])

    # Create an empty raster layer (100 * 100 by default)
    if (verbose) {
      cat("  Use smoothing = \"equal\"\n")
      msg <- paste("  Create an empty raster layer: nrow = ", nrow, 
                   ", ncol = ", ncol, " ...\n", sep = "")
      cat(msg)
      begTimeSub <- Sys.time()
    }
    r1 <- raster(nrows = nrow, ncols = ncol, 
                 xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx)

    # Rasterize the input data set (i.e., polygons) using the empty raster layer
    # you just created. Depending on the input data set, the object 'r2' may
    # contain many empty (NA) cells (e.g., the Auckland metropolitan areas have
    # 5960 empty cells and 4940 non-NA cells). The object 'studyarea' stores the
    # cell IDs with non-NA values. This ensures that we can correctly smooth the 
    # population over the geographic units.
    if (verbose) {
      msg <- paste("  DONE! [", 
                   as.numeric(difftime(Sys.time(), begTimeSub, units = "sec")), 
                   " seconds]\n", sep = "")
      cat(msg)
      cat("  Rasterise 'x' - this may take some time ...\n")
      begTimeSub <- Sys.time()
    }
    r2 <- rasterize(x, r1, field = 0)
    studyarea <- which(!is.na(r2[]))

    # Transform the rasterized surface into a point data set
    if (verbose) {
      msg <- paste("  DONE! [", 
                   as.numeric(difftime(Sys.time(), begTimeSub, units = "sec")), 
                   " seconds]\n", sep = "")
      cat(msg)
      cat("  Transform the rasterised surface into a point data set ...\n")
      begTimeSub <- Sys.time()
    }
    coords <- rasterToPoints(r2)

    # Re-distribute the population counts
    if (verbose) {
      msg <- paste("  DONE! [", 
                   as.numeric(difftime(Sys.time(), begTimeSub, units = "sec")), 
                   " seconds]\n", sep = "")
      cat(msg)
      cat("  Re-distribute the population counts ...\n")
      begTimeSub <- Sys.time()
    }   
    INDEX <- coords[,3]; coords <- coords[,1:2]
    cellPerPolygon <- as.vector(table(INDEX))
    data <- apply(data, 2, function(z) z/cellPerPolygon)
    data <- data[INDEX,]
    # cat(paste("nrow =", nrow(data), "\n"))
    # cat(paste("ncol =", ncol(data), "\n"))

    # Estimate the processing time
    if (verbose) {
      msg <- paste("  DONE! [", 
                   as.numeric(difftime(Sys.time(), begTimeSub, units = "sec")), 
                   " seconds]\n", sep = "")
      cat(msg)
      msg <- paste("DONE! [", 
                   as.numeric(difftime(Sys.time(), begTime, units = "sec")), 
                   " seconds]\n\n", sep = "")
      cat(msg)
    }
    
    for (i in 1:ncol(data)) {
      smoothed.df[[i]] <- r2
      smoothed.df[[i]][studyarea] <- data[,i]
    }
  }
  
  # ----------------------------------------------------------------------------
  # (2) Kernel smoothing
  # ----------------------------------------------------------------------------
  else if (smoothing == "kernel") {
    if (verbose)
      cat("  Use smoothing = \"kernel\"\n")
    
    # If 'window' is not specified, create one for the smoothing procedures
    if (missing(window)) {
      xrange <- range(coords[,1])
      yrange <- range(coords[,2])
      window <- owin(xrange, yrange)
    }
    xPPP <- ppp(coords[,1], coords[,2], window = window, marks = data)

    # Compute a kernel smoothed intensity function. NOTE that in version 0.1,
    # smooth.ppp() was used to construct a smooth surface. It seems, however,
    # that the use of density.ppp() is more appropriate for what I am really 
    # trying to achieve here. Need more thoughts though.
    if (verbose) {
      cat("  Start spatial smoothing of 'x' ...\n")
      begTimeSub <- Sys.time()
    }
    if (missing(sigma))
      sigma <- as.numeric(bw.diggle(xPPP))
    kernl <- list()
    for (i in 1:ncol(xPPP$marks)) {
      kernl[[i]] <- density.ppp(xPPP, w = window, sigma = sigma, 
                                weights = xPPP$marks[,i], dimyx = c(ncol, nrow))
      kernl[[i]][] <- kernl[[i]][] * sum(xPPP$marks[,i])/sum(kernl[[i]][])
    }
    # Transform the rasterized surface into a point data set
    if (verbose) {
      msg <- paste("  DONE! [", 
                   as.numeric(difftime(Sys.time(), begTimeSub, units = "sec")), 
                   " seconds]\n", sep = "")
      cat("  Transform the rasterised surface into a point data set ...\n")
      begTimeSub <- Sys.time()
    }
    for (i in 1:length(kernl)) {
      tmp <- as(kernl[[i]], "RasterLayer")
      if (i == 1) {
        coords <- rasterToPoints(tmp)[,1:2]
        dataNew <- matrix(ncol = length(kernl), nrow = nrow(coords))
      }
      dataNew[,i] <- tmp[]
      # This process may take some time, so let the user know how many columns
      # (i.e., variables) have been processed and update every iteration.
      if (verbose)
        cat("    Column", i, "is processed\n")
    }
    data <- dataNew
    
    # Estimate the processing time
    if (verbose) {
      msg <- paste("  DONE! [", 
                   as.numeric(difftime(Sys.time(), begTimeSub, units = "sec")), 
                   " seconds]\n", sep = "")
      cat(msg)
      msg <- paste("DONE! [", 
                   as.numeric(difftime(Sys.time(), begTime, units = "sec")), 
                   " seconds]\n\n", sep = "")
      cat(msg)
    }
    smoothed.df <- kernl; rm(kernl)
  }
    
  ##############################################################################
  #
  # Step 2. Calculate the population composition of the local environment
  #
  ##############################################################################
  if (verbose) {
    begTime <- Sys.time()
    cat("Estimate the local environment of 'x' - this may take some time ...\n")
  }
  env <- getSegLocalEnv(coords, data, ...)

  if (verbose) {
    msg <- paste("DONE! [", 
                 as.numeric(difftime(Sys.time(), begTime, units = "sec")), 
                 " seconds]\n\n", sep = "")
    cat(msg)
    begTime <- Sys.time()
    cat("Calculate the spatial segregation indices ...\n")
  }

  ##############################################################################
  #
  # Step 3. Compute the segregation indices
  #
  ##############################################################################
  results <- SegSpatial(env, method, useC, negative.rm, tol)

  # VERBOSE --------------------------------------------------------------------
  if (verbose) {
    msg <- paste("DONE! [", 
                 as.numeric(difftime(Sys.time(), begTime, units = "sec")), 
                 " seconds]\n\n", sep = "")
    cat(msg)
    msg <- paste("Total running time: ", 
                 as.numeric(difftime(Sys.time(), excTime, units = "sec")), 
                 " seconds.\n", sep = "")
    cat(msg)
  } # --------------------------------------------------------------------------

  if (missing(sigma))
    sigma <- numeric()
  
  new("SegSpatialExt", nrow = nrow, ncol = ncol, method = method, 
                       sigma = sigma, localenv = env, smooth = smoothed.df,
                       d = results@d, r = results@r, 
                       h = results@h, p = results@p)
}



# ------------------------------------------------------------------------------
# Internal function '.SEGDATA'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
.SEGDATA <- function(x, data, verbose) {
  # VERBOSE --------------------------------------------------------------------
  if (verbose) {
    begTime <- Sys.time()
    cat("Processing 'x' ...\n")
  }

  ##############################################################################
  #
  # (1) If 'x' is an object of class "Spatial" or one that extends the class
  # (e.g., SpatialPoints), then do the following:
  #
  ##############################################################################
  if (inherits(x, "Spatial")) {
    # VERBOSE ------------------------------------------------------------------
    if (verbose)
      cat("  'x' is an object of class \"Spatial\"\n") # -----------------------
    
    coords <- try(coordinates(x), silent = TRUE)
    if (class(coords) == "try-error")
      stop("failed to extract coordinates from 'x'", call. = FALSE)

    # VERBOSE ------------------------------------------------------------------
    if (verbose) {
      msg <- paste("  ", nrow(coords), " coordinates extracted", sep = "")
      msg <- paste(msg, " from 'x' succesfully\n", sep = "")
      cat(msg)
    } # ------------------------------------------------------------------------
    
    if (missing(data)) {
      # The code below will success if 'x' includes a data frame (e.g., 'x' is
      # a SpatialPointsDataFrame object).
      data <- try(as.matrix(x@data), silent = TRUE) 
      if (class(data) == "try-error")
        stop("'data' is missing, with no default", call. = FALSE)

      # VERBOSE ----------------------------------------------------------------
      if (verbose)
        cat("  'data' is missing, use the one attached to 'x'\n") # ------------
    } else {
      data <- as.matrix(data)
    }
    
    # VERBOSE ------------------------------------------------------------------
    if (verbose)
      cat("  check if 'data' has any NA values\n") # ---------------------------
    removeNA <- which(apply(data, 1, function(z) any(is.na(z))))
    if (length(removeNA) > 0) {
      data <- data[-removeNA,]
      coords <- coords[-removeNA,]
      # VERBOSE ----------------------------------------------------------------
      if (verbose) {
        msg <- paste(length(removeNA), "NA(s) found and removed\n")
        cat(msg)
      } # ----------------------------------------------------------------------
    }
  }

  ##############################################################################
  #
  # (2) If 'x' is an object of class "ppp", then do the following:
  #
  ##############################################################################
  else if (is(x, "ppp")) {
    # VERBOSE ------------------------------------------------------------------
    if (verbose)
      cat("  'x' is an object of class \"ppp\"\n") # ---------------------------
    
    coords <- cbind(x = x$x, y = x$y)
    if (missing(data)) {
      # The code below will success if 'x' includes a data frame (i.e., marks).
      if (is.null(x$marks))
        stop("'data' is missing, with no default", call. = FALSE)
      else {
        data <- as.matrix(x$marks)
        # VERBOSE --------------------------------------------------------------
        if (verbose)
          cat("  'data' is missing, use the one attached to 'x'\n") # ----------
      }
    } else {
      data <- as.matrix(data)
    }
  }

  ##############################################################################
  #
  # (3) If 'x' is a n * 2 matrix or data frame object, then do the following:
  #
  ##############################################################################
  else if (is.matrix(x) || is.data.frame(x)) {
    # VERBOSE ------------------------------------------------------------------
    if (verbose)
      cat("  'x' is an object of class \"matrix\" or \"data.frame\"\n") # ------
    
    coords <- as.matrix(x)
    if (ncol(coords) != 2 || !is.numeric(coords))
      stop("'x' must be a numeric matrix with two columns", call. = FALSE)
    if (missing(data))
      stop("'data' is missing, with no default", call. = FALSE)
    else
      data <- as.matrix(data)
  }

  ##############################################################################
  #
  # (4) If 'x' is not one of the supporting classes:
  #
  ##############################################################################
  else
    stop("invalid object 'x'", call. = FALSE)

  # VERBOSE --------------------------------------------------------------------
  if (verbose) {
    msg <- paste("DONE! [", 
                 as.numeric(difftime(Sys.time(), begTime, units = "sec")), 
                 " seconds]\n\n", sep = "")
    cat(msg)
  }
  
  invisible(list(coords = coords, data = data))
}
