\name{seg}
\alias{seg}
\title{Calculate Index of Dissimilarity}
\description{
\code{seg} calculates Duncan and Duncan's index of dissimilarity between two population groups. If \sQuote{nb} is given, the index can be adjusted to reflect the spatial distribution of population.
}
\usage{
seg(data, nb, verbose = FALSE)
}
\arguments{
  \item{data}{a numeric matrix or data frame with two columns that represent mutually exclusive population groups (e.g., Asians and non-Asians). If more than two columns are given, only the first two will be used for computing the index.}
  \item{nb}{an optional \code{matrix} object describing the intensity of interaction between geographic units.}
  \item{verbose}{logical. If TRUE, print the current stage of the computation and time spent on each job to the screen.}
}
\details{
\sQuote{nb} must be a square matrix (same number of rows and columns) but does not have to be symmetric. When \sQuote{nb} is not specified, \code{seg} calculates the traditional index of dissimilarity proposed by Duncan and Duncan (1955). 

If \sQuote{nb} is a rook-based contiguity matrix standardised by the total number of neighbours, the calculation follows the Morrill's adjusted index of dissimilarity (1991). See the example code below and \code{help(nb2mat)} for more information regarding how to construct such a matrix. 

If \sQuote{nb} contains row-standardised shared boundary lengths between geographic units, the computed index is Wong's adjusted D. In \R, \code{vect2neigh} in the \pkg{spgrass6} package can be used to obtain the shared boundary lengths.
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
\code{\link{spseg}}, \code{\link{whiteseg}}
}
\examples{
if (require(spdep)) {
  data(segdata)
  segdata.df <- data.frame(segdata)
  segdata.nb <- nb2mat(poly2nb(segdata, queen = FALSE), style = "B")
  segdata.nb <- segdata.nb / sum(segdata.nb)
  dd <- numeric(); dm <- numeric()
  pd <- par()
  par(mfrow = c(2, 4), mar = c(0, 1, 0, 1))
  for (i in 1:8) {
    idx <- 2 * i
    dd <- append(dd, seg(segdata.df[,(idx-1):idx]))
    dm <- append(dm, seg(segdata.df[,(idx-1):idx], segdata.nb))
    full <- segdata.df[,(idx-1)] == 100
    half <- segdata.df[,(idx-1)] == 50
    plot(segdata)
    plot(segdata[full,], col = "Black", add = TRUE)
    if (any(half))
      plot(segdata[half,], col = "Grey", add = TRUE)
    text(5, 11.5, labels = paste("D = ", round(dd[i], 2), 
                                 ", D(adj) = ", round(dm[i], 2), sep = ""))
  }
}
}
