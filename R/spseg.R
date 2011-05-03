# ------------------------------------------------------------------------------
# Function 'spseg'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
spseg <- 
  function(x, data, method = "all", smoothing = "none", nrow = 100, ncol = 100, 
           window, sigma, useC = TRUE, verbose = FALSE, ...) {
  
  ### verbose = TRUE ###########################################################
  if (verbose) {
    excTime <- Sys.time()
    begTime <- Sys.time()
    cat("Processing 'x' ...\n")
  }
  ##############################################################################
  
  # If 'x' is an object of class "Spatial" or one that extends the class
  # (e.g., SpatialPoints), then do the following:
  if (inherits(x, "Spatial")) {
  
    ### verbose = TRUE #########################################################
    if (verbose) {
      cat("  'x' is an object of class \"Spatial\"\n")
    }
    ############################################################################
    
    coords <- try(coordinates(x), silent = TRUE)
    if (class(coords) == "try-error")
      stop("failed to extract coordinates from 'x'", call. = FALSE)

    ### verbose = TRUE #########################################################
    if (verbose) {
      msg <- paste("  ", nrow(coords), " coordinates extracted", sep = "")
      msg <- paste(msg, " from 'x' succesfully\n", sep = "")
      cat(msg)
    }
    ############################################################################
    
    if (missing(data)) {
      # The code below will success if 'x' includes a data frame (e.g., 'x' is
      # a SpatialPointsDataFrame object).
      data <- try(as.matrix(x@data), silent = TRUE) 
      if (class(data) == "try-error")
        stop("'data' is missing, with no default", call. = FALSE)

      ### verbose = TRUE #######################################################
      if (verbose) {
        cat("  'data' is missing, use the one attached to 'x'\n")
      }
      ##########################################################################
      
    } else {
      data <- as.matrix(data)
      
      ### verbose = TRUE #######################################################   
      if (verbose) {
        cat("  make sure 'data' is an object of class \"matrix\"\n")
      }
      ##########################################################################
    }
    
    ### verbose = TRUE #########################################################
    if (verbose) {
      cat("  check if 'data' has any NA values\n")
    }
    ############################################################################
    
    removeNA <- which(apply(data, 1, function(z) any(is.na(z))))
    if (length(removeNA) > 0) {
      data <- data[-removeNA,]
      coords <- coords[-removeNA,]
      
      ### verbose = TRUE #######################################################
      if (verbose) {
        msg <- paste(length(removeNA), "NA found and removed\n")
        cat(msg)
      }
      ##########################################################################
    }
  } 

  # If 'x' is an object of class "ppp", then do the following:
  else if (is(x, "ppp")) {
  
    ### verbose = TRUE #########################################################
    if (verbose) {
      cat("  'x' is an object of class \"ppp\"\n")
    }
    ############################################################################
    
    coords <- cbind(x = x$x, y = x$y)
    if (missing(data)) {
      if (x$mark != "none")
        stop("'data' is missing, with no default", call. = FALSE)
      else
        data <- as.matrix(x$mark)
        
      ### verbose = TRUE #######################################################
      if (verbose) {
        cat("  'data' is missing, use the one attached to 'x'\n")
      }
      ##########################################################################
    } else {
      data <- as.matrix(data)
      
      ### verbose = TRUE #######################################################
      if (verbose) {
        cat("  make sure 'data' is an object of class \"matrix\"\n")
      }
      ##########################################################################
    }
  }
  
  # If 'x' is a n * 2 matrix or data frame object, then do the following:
  else if (is.matrix(x) || is.data.frame(x)) {
  
    ### verbose = TRUE #########################################################
    if (verbose) {
      cat("  'x' is an object of class \"matrix\" or \"data.frame\"\n")
    }
    ############################################################################
    
    coords <- as.matrix(x)
    if (ncol(coords) != 2 || !is.numeric(coords))
      stop("'x' must be a numeric matrix with two columns", call. = FALSE)
    if (missing(data))
      stop("'data' is missing, with no default", call. = FALSE)
    else
      data <- as.matrix(data)
      
    ### verbose = TRUE #########################################################
    if (verbose) {
      cat("  make sure 'data' is an object of class \"matrix\"\n")
    }
    ############################################################################
  }
  
  # If 'x' is not one of the supporting classes:
  else {
    stop("invalid object 'x'", call. = FALSE)
  }


  ### verbose = TRUE ###########################################################
  if (verbose) {  ## Verbose
    msg <- paste("DONE! [", 
                 as.numeric(difftime(Sys.time(), begTime, units = "sec")), 
                 " seconds]\n\n", sep = "")
    cat(msg)
  }
  ##############################################################################


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

  # ----------------------------------------------------------------------------
  #
  # Step 1. Estimate the data surface
  #
  # ----------------------------------------------------------------------------  

  ### verbose = TRUE ###########################################################
  if (verbose) {
    begTime <- Sys.time()
    cat("Create a data surface of 'x' ...\n")
  }
  ##############################################################################

  # (1) If 'smoothing' is "equal":
  if (smoothing == "equal") {
    if (!inherits(x, "SpatialPolygons")) {
      msg <- paste("'x' should be a SpatialPolygons object to", 
                   "use the \"equal\" smoothing option")
      stop(msg, call. = FALSE)
    }

    xmn <- min(coords[,1]); xmx <- max(coords[,1])
    ymn <- min(coords[,2]); ymx <- max(coords[,2])

    ### verbose = TRUE #########################################################
    if (verbose) {
      cat("  Use smoothing = \"equal\"\n")
      msg <- paste("  Create an empty raster layer: nrow = ", nrow, 
                   ", ncol = ", ncol, " ...\n", sep = "")
      cat(msg)
      begTimeSub <- Sys.time()
    }
    ############################################################################
    
    r1 <- raster(nrows = nrow, ncols = ncol, 
                 xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx)

    ### verbose = TRUE #########################################################
    if (verbose) {
      msg <- paste("  DONE! [", 
                   as.numeric(difftime(Sys.time(), begTimeSub, units = "sec")), 
                   " seconds]\n", sep = "")
      cat(msg)
      cat("  Rasterise 'x' - this may take some time ...\n")
      begTimeSub <- Sys.time()
    }
    ############################################################################
    
    r2 <- rasterize(x, r1, field = 0)

    ### verbose = TRUE #########################################################
    if (verbose) {
      msg <- paste("  DONE! [", 
                   as.numeric(difftime(Sys.time(), begTimeSub, units = "sec")), 
                   " seconds]\n", sep = "")
      cat(msg)
      cat("  Transform the rasterised surface into a point data set ...\n")
      begTimeSub <- Sys.time()
    }
    ############################################################################
    
    coords <- rasterToPoints(r2)

    ### verbose = TRUE #########################################################
    if (verbose) {
      msg <- paste("  DONE! [", 
                   as.numeric(difftime(Sys.time(), begTimeSub, units = "sec")), 
                   " seconds]\n", sep = "")
      cat(msg)
      cat("  Re-distribute the population counts ...\n")
      begTimeSub <- Sys.time()
    }
    ############################################################################
    
    INDEX <- coords[,3]; coords <- coords[,1:2]
    cellPerPolygon <- as.vector(table(INDEX))
    data <- apply(data, 2, function(z) z/cellPerPolygon)
    data <- data[INDEX,]

    ### verbose = TRUE #########################################################
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
    ############################################################################
  }
  
  # (2) If 'smoothing' is "kernel":
  else if (smoothing == "kernel") {
  
    ### verbose = TRUE #########################################################
    if (verbose) {
      cat("  Use smoothing = \"kernel\"\n")
    }
    ############################################################################
    
    if (missing(window)) {
      xrange <- range(coords[,1])
      yrange <- range(coords[,2])
      window <- owin(xrange, yrange)
    }

    xPPP <- ppp(coords[,1], coords[,2], window = window, marks = data)

    ### verbose = TRUE #########################################################
    if (verbose) {
      cat("  Start spatial smoothing of 'x' ...\n")
      begTimeSub <- Sys.time()
    }
    ############################################################################

    if (missing(sigma))      
      kernl <- smooth.ppp(xPPP, w = window, dimyx = c(ncol, nrow))
    else
      kernl <- smooth.ppp(xPPP, w = window, sigma, dimyx = c(ncol, nrow))

    ### verbose = TRUE #########################################################
    if (verbose) {
      msg <- paste("  DONE! [", 
                   as.numeric(difftime(Sys.time(), begTimeSub, units = "sec")), 
                   " seconds]\n", sep = "")
      cat("  Transform the rasterised surface into a point data set ...\n")
      begTimeSub <- Sys.time()
    }
    ############################################################################

    for (i in 1:length(kernl)) {
      tmp <- as(kernl[[i]], "RasterLayer")
      if (i == 1) {
        coords <- rasterToPoints(tmp)[,1:2]
        dataNew <- matrix(ncol = length(kernl), nrow = nrow(coords))
      }
      dataNew[,i] <- tmp[]

      ### verbose = TRUE #######################################################
      if (verbose) {
        cat("    Column", i, "is processed\n")
      }
      ##########################################################################
      
    }
    data <- dataNew
    
    ### verbose = TRUE #########################################################
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
    ############################################################################
    
  }
    
  # ----------------------------------------------------------------------------
  #
  # Step 2. Calculate the population composition of the local environment
  #
  # ----------------------------------------------------------------------------

  ### verbose = TRUE ###########################################################
  if (verbose) {
    begTime <- Sys.time()
    cat("Estimate the local environment of 'x' - this may take some time ...\n")
  }
  ##############################################################################

  env <- getSegLocalEnv(coords, data, ...)


  ### verbose = TRUE ###########################################################
  if (verbose) {
    msg <- paste("DONE! [", 
                 as.numeric(difftime(Sys.time(), begTime, units = "sec")), 
                 " seconds]\n\n", sep = "")
    cat(msg)
    begTime <- Sys.time()
    cat("Calculate the spatial segregation indices ...\n")
  }
  ##############################################################################

  # ----------------------------------------------------------------------------
  #
  # Step 3. Compute the segregation indices
  #
  # ----------------------------------------------------------------------------
  results <- SegSpatial(env, method, useC)

  ### verbose = TRUE ###########################################################
  if (verbose) {
    msg <- paste("DONE! [", 
                 as.numeric(difftime(Sys.time(), begTime, units = "sec")), 
                 " seconds]\n\n", sep = "")
    cat(msg)
    msg <- paste("Total running time: ", 
                 as.numeric(difftime(Sys.time(), excTime, units = "sec")), 
                 " seconds.\n", sep = "")
    cat(msg)
  }
  ##############################################################################

  results
}

