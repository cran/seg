\name{spseg}
\alias{spseg}
\alias{SegSpatial}
\title{Calculate Spatial Segregation Measures}
\description{
Calculate Reardon and O'Sullivan's four spatial segregation measures.
}
\usage{
spseg(x, data, method = "all", smoothing = "none", nrow = 100, 
      ncol = 100, window, sigma, useC = TRUE, verbose = FALSE, \dots)
SegSpatial(env, method = "all", useC = TRUE)
}
\arguments{
  \item{x}{a numeric matrix or data frame with coordinates (each row is a point), or an object of class \code{Spatial} or \code{ppp}.}
  \item{env}{an object of class \code{SegLocalEnv}.}
  \item{data}{an object of class \code{matrix}, or one that can be coerced to that class. The number of rows in \sQuote{data} should equal the number of points in \sQuote{x}, and the number of columns should be greater than one (at least two population groups are required). This can be missing if \sQuote{x} has a data frame attached to it.}
  \item{method}{a vector of character strings indicating an measure or measures to be computed. This must be one or more of the strings "all" (default), "exposure", "information", "diversity", and "dissimilarity". Abbreviations are accepted, as long as it is clear which method is meant.}
  \item{smoothing}{a character string indicating a smoothing method. This must be (an abbreviation of) one of the strings "none" (default), "kernel", or "equal".}
  \item{nrow}{an optional numeric value indicating the number of row cells in the rasterised data surface. Ignored if \sQuote{smoothing = "none"}.}
  \item{ncol}{an optional numeric value indicating the number of column cells.}
  \item{window}{an optional object of class \code{owin} to be passed to \code{smooth.ppp}. See \sQuote{Details}.}
  \item{sigma}{an optional numeric value specifying standard deviation of isotropic Gaussian smoothing kernel to be passed to \code{smooth.ppp}. See also \sQuote{Details}.}
  \item{useC}{logical. If TRUE, calculate the segregation values in C.}
  \item{verbose}{logical. If TRUE, print the current stage of the computation and time spent on each job to the screen.}
  \item{\dots}{optional arguments to be passed to \code{getLocalEnv} to compute the population composition of each local environment.}
}
\details{
\code{SegSpatial} computes the set of spatial segregation measures proposed by Reardon and O'Sullivan.

\code{spseg} is a wrapper function, which calls \code{SegSpatial} after constructing a population density surface and its local environment parameters with user-specified options. Currently the population density surface is estimated using the \code{rasterize} function in the \code{raster} package if the population density is assumed to be uniform in each census tract, or using \code{smooth.ppp} in the \code{spatstat} package if the kernel density estimation is to be used. The local environment parameters are calculated based on the output surface using \code{\link{getSegLocalEnv}}.
}
\value{
An object of class \code{SegSpatial}.
}
\note{
The exposure/isolation index, P, is presented in a matrix form. The spatial exposure of group \sQuote{m} to group \sQuote{n} is located in the row \sQuote{m} and column \sQuote{n} of the matrix. The matrix is rarely symmetric in practice so the spatial exposure index should be interpreted with care.

The spatial isolation index values are given in the diagonal cells of the matrix; cell value at (m, m) indicates the degree of spatial isolation for group \sQuote{m} for example.
}
\author{
Seong-Yun Hong
}
\references{
Reardon, S. F. and O'Sullivan, D. (2004) Measures of spatial segregation. \emph{Sociological Methodology}, \bold{34}, 121-162.

Reardon, S. F., Farrell, C. R., Matthews, S. A., O'Sullivan, D., Bischoff, K., and Firebaugh, G. (2009) Race and space in the 1990s: Changes in the geographic scale of racial residential segregation, 1990-2000. \emph{Social Science Research}, \bold{38}, 55-70.
}
\seealso{
\code{\link{getSegLocalEnv}}, \code{\link{SegSpatial-class}}, \code{rasterize}, \code{smooth.ppp}
}
\examples{
# Create a random data set with 50 data points and 3 population groups
xy <- matrix(runif(100), ncol = 2)
pop <- matrix(runif(150), ncol = 3)
r1 <- spseg(xy, pop, smoothing = "kernel", maxdist = 0.5)
r2 <- spseg(xy, pop, smoothing = "kernel", useExp = FALSE, 
            power = 0, maxdist = 0.5)
print(r2, digits = 3)

# Auckland population data set
data(auckpop2006)
r3 <- spseg(auckpop2006, smoothing = "kernel", maxdist = 3000)
\dontrun{
r4 <- spseg(auckpop2006, smoothing = "equal", maxdist = 100)
r5 <- spseg(auckpop2006, smoothing = "equal", maxdist = 3000)
r6 <- spseg(auckpop2006, smoothing = "kernel", nrow = 300, ncol = 300, 
            maxdist = 3000, verbose = TRUE)
}
}
