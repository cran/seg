\name{whiteseg}
\alias{whiteseg}
\title{Calculate White's Segregation Measure}
\description{
\code{whiteseg} computes White's P index. Unlike the traditional, aspatial index of dissimilarity, this measure incorporates the spatial relationships among the geographic units into the calculation, so it does not suffer from the checkerboard problem.
}
\usage{
whiteseg(x, data, nb, fun, ...)
}
\arguments{
  \item{x}{a numeric matrix or data frame with coordinates (each row is a point), or an object of class \code{Spatial}.}
  \item{data}{an object of class \code{matrix}, or one that can be coerced to that class. The number of rows in \sQuote{data} should equal the number of geographic units in \sQuote{x}, and the number of columns should be greater than one (at least two population groups are required). This can be missing if \sQuote{x} has a data frame attached to it.}
  \item{nb}{an optional \code{matrix} object indicating the distances between the geographic units.}
  \item{fun}{a function for the calculation of proximity. The function should take a numeric vector as an argument (distance) and return a vector of the same length (proximity). If this is not specified, a negative exponential function is used by default.}
  \item{...}{optional arguments to be passed to \code{dist} when calculating the distances between the geographic units in \sQuote{x}. Ignored if \sQuote{nb} is given. See \code{help(dist)} for available options.}
}
\details{
\sQuote{nb} must be a square matrix (same number of rows and columns) but does not have to be symmetric. When \sQuote{nb} is not specified, \code{whiteseg} computes the distance matrix of \sQuote{x} using the specified optional arguments in \sQuote{...}.
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
require(sp)
grd <- GridTopology(cellcentre.offset=c(0.5,0.5),
                    cellsize=c(1,1), cells.dim=c(10,10))
grd.sp <- as.SpatialPolygons.GridTopology(grd)

data(segdata)
parval <- par()
par(mfrow = c(2, 4), mar = c(0, 1, 0, 1))

d <- numeric()
p <- numeric()
for (i in 1:8) {
  idx <- 2 * i
  d <- append(d, seg(segdata[,(idx-1):idx]))
  p <- append(p, whiteseg(grd.sp, data = segdata[,(idx-1):idx]))
  full <- segdata[,(idx-1)] == 100
  half <- segdata[,(idx-1)] == 50
  plot(grd.sp)
  plot(grd.sp[full,], col = "Black", add = TRUE)
  if (any(half))
    plot(grd.sp[half,], col = "Grey", add = TRUE)
  text(5, 11.5, labels = 
       paste("D = ", round(d[i], 2), ", P = ", round(p[i], 2), sep = ""))
}
par(mfrow = parval$mfrow, mar = parval$mar)
}
