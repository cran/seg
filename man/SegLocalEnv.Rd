\name{SegLocalEnv}
\alias{SegLocalEnv}
\title{Create an Object of Class SegLocalEnv}
\description{
Create a new object of class \code{SegLocalEnv} from a matrix of coordinates and population data.
}
\usage{
SegLocalEnv(coords, data, env, proj4string = CRS(as.character(NA)))
}
\arguments{
  \item{coords}{a numeric matrix or data frame with coordinates (each row is a point).}
  \item{data}{an object of class \code{matrix} containing the population of each group at each data point. The number of rows in \sQuote{data} should equal the number of points in \sQuote{coords}, and the number of columns should be greater than one (at least two population groups are required).}
  \item{env}{an optional object of class \code{matrix} containing the local environment parameters. Must be the same dimensions as \sQuote{data}. If missing, use \sQuote{data}.}
  \item{proj4string}{an optional projection string of class \code{CRS}.}
}
\value{
An object of class \code{SegLocalEnv}.
}
\author{
Seong-Yun Hong
}
\seealso{
\code{\link{getSegLocalEnv}}, \code{\link{update.SegLocalEnv}}, \code{\link{spseg}}, \code{\link{SegLocalEnv-class}}
}
\examples{
# Create a random data set with 50 data points and 3 population groups
xy <- matrix(runif(100), ncol = 2)
pop <- matrix(runif(150), ncol = 3)
localenv <- SegLocalEnv(coords = xy, data = pop)
}
