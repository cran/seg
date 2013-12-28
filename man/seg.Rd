\name{seg}
\alias{seg}
\title{Index of Dissimilarity}
\description{
Calculates the index of dissimilarity proposed by Duncan and Duncan (1955). If \sQuote{nb} is given, the index is adjusted to reflect the spatial distribution of population.
}
\usage{
seg(data, nb)
}
\arguments{
  \item{data}{a numeric matrix or data frame with two columns that represent mutually exclusive population groups (e.g., Asians and non-Asians). If more than two columns are given, only the first two will be used for computing the index.}
  \item{nb}{an optional \code{matrix} object describing the intensity of interaction between geographic units.}
}
\details{
\sQuote{nb} must be a square matrix (i.e., one that has the same number of rows and columns) but does not have to be symmetric. If \sQuote{nb} is not specified, \code{seg} calculates the traditional index of dissimilarity proposed by Duncan and Duncan (1955). 

If \sQuote{nb} is a rook-based contiguity matrix standardised by the total number of neighbours, the index is adjusted as described in Morrill (1991). See the example code below and \code{nb2mat} in \pkg{spdep} for more information regarding how to construct such a matrix. 

For Wong's D(w) and D(s), see \url{https://sites.google.com/site/hongseongyun/seg}. Note that the sum of all elements in \sQuote{nb} should equal one.
}
\value{
A single numeric value between 0 and 1, indicating the degree of segregation; 0 for no segregation, and 1 for complete segregation.
}
\author{
Seong-Yun Hong
}
\references{
Duncan, O. D., & Duncan, B. (1955). A methodological analysis of segregation indexes. \emph{American Sociological Review}, \bold{20}, 210-217.

Morrill, R. L. (1991). On the measure of geographic segregation. \emph{Geography Research Forum}, \bold{11}, 25-36.

Wong, D. W. S. (1993). Spatial indices of segregation. \emph{Urban Studies}, \bold{30}, 559-572.
}
\seealso{
\code{\link{spseg}}, \code{nb2mat}
}
\examples{
if (require(spdep)) { # package 'spdep' is required

  # uses the idealised landscapes in 'segdata'
  data(segdata)
  grd <- GridTopology(cellcentre.offset=c(0.5,0.5),
                      cellsize=c(1,1), cells.dim=c(10,10))
  grd.sp <- as.SpatialPolygons.GridTopology(grd)
  grd.nb <- nb2mat(poly2nb(grd.sp, queen = FALSE), style = "B")
  grd.nb <- grd.nb / sum(grd.nb)

  d <- numeric(); m <- numeric()
  for (i in 1:8) {
    idx <- 2 * i
    d <- append(d, seg(segdata[,(idx-1):idx]))
    m <- append(m, seg(segdata[,(idx-1):idx], grd.nb))
    full <- segdata[,(idx-1)] == 100
    half <- segdata[,(idx-1)] == 50
    plot(grd.sp)
    plot(grd.sp[full,], col = "Black", add = TRUE)
    if (any(half))
      plot(grd.sp[half,], col = "Grey", add = TRUE)
    text(5, 11.5, labels = paste("D = ", round(d[i], 2), 
                                 ", D(adj) = ", round(m[i], 2), sep = ""))
  }
}
}
