\name{whiteseg}
\alias{whiteseg}
\title{Calculate White's Segregation Measure}
\description{
\code{whiteseg} computes White's segregation statistic. Unlike the traditional index of dissimilarity, this is a spatial measure that takes the spatial arrangements of population into account.
}
\usage{
whiteseg(x, data, nb, method = "euclidean", p = 2, fun, verbose = FALSE)
}
\arguments{
  \item{x}{a numeric matrix or data frame with coordinates (each row is a point), or an object of class \code{Spatial}.}
  \item{data}{an object of class \code{matrix}, or one that can be coerced to that class. The number of rows in \sQuote{data} should equal the number of geographic units in \sQuote{x}, and the number of columns should be greater than one (at least two population groups are required). This can be missing if \sQuote{x} has a data frame attached to it.}
  \item{nb}{an optional \code{matrix} object indicating the distances between the geographic units.}
  \item{method}{a length-one character string indicating the distance measure to be used. This argument is passed to \code{dist} to calculate the distances between the geographic units in \sQuote{x}. See \code{help(dist)} for available options.}
  \item{p}{an optional argument to be passed to \code{dist}.}
  \item{fun}{a function for the calculation of proximity. The function should take a numeric vector as an argument (distance) and return a vector of the same length (proximity). If this is not specified, a negative exponential function is used by default.}
  \item{verbose}{logical. If TRUE, print the current stage of the computation and time spent on each job to the screen.}
}
\details{
\sQuote{nb} must be a square matrix (same number of rows and columns) but does not have to be symmetric. When \sQuote{nb} is not specified, \code{whiteseg} computes the distance matrix of \sQuote{x} using the specified optional arguments \sQuote{method} and \sQuote{p}.
}
\value{
A single numeric value indicating the degree of segregation; a value of 1 indicates absence of segregation, and values greater than 1.0 indicate clustering. If the index value is less than one, it indicates an unusual form of segregation (i.e., people live closer to other population groups).
}
\author{
Seong-Yun Hong
}
\references{
White, M. J. (1983). The measurement of spatial segregation. \emph{The American Journal of Sociology}, \bold{88}, 1008-1018.
}
\seealso{
\code{\link{seg}}, \code{\link{spseg}}, \code{dist}
}
\examples{
if (require(spdep)) {
  data(segdata)
  segdata.df <- data.frame(segdata)
  segdata.nb <- nb2mat(poly2nb(segdata, queen = FALSE), style = "B")
  segdata.nb <- segdata.nb / sum(segdata.nb)
  d <- numeric(); p <- numeric()
  par(mfrow = c(2, 4), mar = c(0, 1, 0, 1))
  for (i in 1:8) {
    idx <- 2 * i
    d <- append(d, seg(segdata.df[,(idx-1):idx], segdata.nb))
    p <- append(p, whiteseg(segdata, data = segdata.df[,(idx-1):idx]))
    full <- segdata.df[,(idx-1)] == 100
    half <- segdata.df[,(idx-1)] == 50
    plot(segdata)
    plot(segdata[full,], col = "Black", add = TRUE)
    if (any(half))
      plot(segdata[half,], col = "Grey", add = TRUE)
    text(5, 11.5, labels = 
      paste("D = ", round(d[i], 2), ", P = ", round(p[i], 2), sep = ""))
  }
  
  data(nzcity)
  dd <- numeric(); pp <- numeric()
  for (i in 1:16) {
    city.df <- data.frame(nzcity[[i]])
    city.df <- cbind(city.df[,4], apply(city.df[,c(2:3, 5:7)], 1, sum))
    city.nb <- nb2mat(poly2nb(nzcity[[i]], queen = FALSE), style = "B",
                      zero.policy = TRUE)
    city.nb <- city.nb / sum(city.nb)
    dd <- append(dd, seg(city.df, city.nb))
    pp <- append(pp, whiteseg(nzcity[[i]], city.df, fun = function(z) 1/z))
  }
  nzcity.lm <- lm(pp ~ dd)
  summary(nzcity.lm)
  plot(dd, pp, xlab = "D(adj)", ylab = "P", 
       main = paste("Residential Segregation of\n",
                    "Pacific Peoples in 16 New Zealand Cities", sep = ""))
  abline(a = nzcity.lm$coefficients[1], 
         b = nzcity.lm$coefficients[2], col = "Red", lty = "dotted")
  which.pos <- c(3, 2, 2, 3, 3, 3, 2, 3, 3, 4, 3, 2, 3, 3, 3, 1)
  text(dd, pp, pos = which.pos, labels = names(nzcity), cex = 0.8)
}
}
