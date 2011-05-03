# ------------------------------------------------------------------------------
# Function 'SegSpatial'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
SegSpatial <- function(env, method = "all", useC = TRUE) {

  validObject(env)
  method <- match.arg(method, c("exposure", "information", "diversity", 
                                "dissimilarity", "all"), several.ok = TRUE)
  if ("all" %in% method)
    method <- c("exposure", "information", "diversity", "dissimilarity")

  P <- matrix(0, nrow = 0, ncol = 0)
  H <- numeric(); R <- numeric(); D <- numeric()
  if (any(env@env == 0))
    warning("zero values will be ignored during the calulcation of H", 
            call. = FALSE)
  
  if (useC) {
    m <- ncol(env@data)
    method <- c("exposure" %in% method, "information" %in% method,
                "diversity" %in% method, "dissimilarity" %in% method)
    tmp <- .Call("spseg", as.vector(env@data), as.vector(env@env), 
                          as.integer(m), as.integer(method))
    results <- list(); n <- m^2
    if (!is.na(tmp[1])) {
      results$p <- matrix(tmp[1:n], ncol = m, byrow = TRUE)
      rownames(results$p) <- colnames(results$p) <- colnames(env@data)
    }
    if (!is.na(tmp[n+1]))
      results$h <- tmp[n+1] 
    if (!is.na(tmp[n+2]))
      results$r <- tmp[n+2]
    if (!is.na(tmp[n+3]))
      results$d <- tmp[n+3]
  } 
  
  else {
    # Number of population groups
    m <- ncol(env@data)
    # Total population in the study area
    ptsSum <- sum(env@data)
    # Population of all groups at each data point
    ptsRowSum <- apply(env@data, 1, sum)
    # Total population of each subgroup
    ptsColSum <- apply(env@data, 2, sum)
    # Proportion of each subgroup
    ptsProp <- ptsColSum / ptsSum
    # Population proportion of each group at each local environment
    envProp <- t(apply(env@env, 1, function(z) z/sum(z)))

    if ("exposure" %in% method) {
      P <- matrix(0, nrow = m, ncol = m)
      rownames(P) <- colnames(P) <- colnames(env@data)
      for (i in 1:m) {
        A <- env@data[, i] / ptsColSum[i]
        for (j in 1:m) {
          P[i, j] <- sum(A * envProp[, j])
        }
      }
    }

    if ("information" %in% method) {
      Ep <- apply(envProp, 1, function(z) sum(z * log(z, base = m))) * -1
      E <- sum(ptsProp * log(ptsProp, base = m)) * -1
      H <- 1 - (sum(ptsRowSum * Ep) / (ptsSum * E))
    }

    if ("diversity" %in% method) {
      Ip <- apply(envProp, 1, function(z) sum(z * (1 - z)))
      I <- sum(ptsProp * (1 - ptsProp))
      R <- 1 - sum((ptsRowSum * Ip) / (ptsSum * I))
    }

    if ("dissimilarity" %in% method) {
      I <- sum(ptsProp * (1 - ptsProp))
      constant <- ptsRowSum / (2 * ptsSum * I)
      Dp <- t(apply(envProp, 1, function(z) abs(z - ptsProp)))
      D <- sum(apply(Dp, 2, function(z) sum(z * constant)))
    }
    results <- list(p = P, h = H, r = R, d = D)
  }
  
  as(results, "SegSpatial")
}
