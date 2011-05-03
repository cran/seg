\name{getSegLocalEnv}
\alias{getSegLocalEnv}
\title{Calculate Local Environment Parameters}
\description{
Calculate the population composition of the local environment at each data point from a matrix of coordinates, or an object of class \code{Spatial} or \code{ppp}.
}
\usage{
getSegLocalEnv(x, data, power = 2, useExp = TRUE, maxdist, 
               sprel, error = .Machine$double.eps)
}
\arguments{
  \item{x}{a numeric matrix or data frame with coordinates (each row is a point), or an object of class \code{Spatial} or \code{ppp}.}
  \item{data}{an object of class \code{matrix}, or one that can be coerced to that class. The number of rows in \sQuote{data} should equal the number of points in \sQuote{x}, and the number of columns should be greater than one (at least two population groups are required). This can be missing if \sQuote{x} has a data frame attached to it.}
  \item{power}{a numeric value that determines the change rate of a distance weight function. If zero, all data points have the same weight regardless of the distance. Typically 1-5.}
  \item{useExp}{logical. If FALSE, use a simple inverse distance function instead of a negative exponential function. See \sQuote{Details}.}
  \item{maxdist}{an optional numeric value specifying a search radius for the construction of each local environment. Must be a positive value, or zero.}
  \item{sprel}{an optional object of class \code{dist} or \code{nb}. See \sQuote{Details}.}
  \item{error}{a numeric value. If \sQuote{useExp = FALSE}, this value will be added to the denominator to prevent dividing by zero.}
}
\details{
The local environment parameters are the weighted averages of the populations of all points within a search radius \sQuote{maxdist} and are an essential component for calculation of the spatial segregation measures.

By default, the weight of each point is calculated from a negative exponential function, which is defined as:

\deqn{w(d) = e^{-d \times power}}{%
      w(d) = exp(-d * power)}

where \emph{d} is the Euclidean distance between two points.

If \sQuote{useExp} is FALSE, a simple inverse distance function is used to calculate the weight of each point:

\deqn{w(d) = \frac{1}{(d + error)^{power}}}{%
      w(d) = 1 / (d + error)^power}

If \sQuote{maxdist} is not provided (default), all data points in the study area are used for the construction of each local environment. It is recommended to specify this parameter to speed up the calculation process.

If a distance measure other than the Euclidean distance is required to represent spatial proximity between the points, the users can provide an object of class \code{dist}, which contains the distances between all pairs of the points, through an optional argument \sQuote{sprel}. One convenient way of obtaining such information may be the use of the function \code{dist}, which offers a variety of distance measures, such as Manhattan, Canberra, and Minkowski.

Or alternatively, one can supply an object of class \code{nb} to use a k-nearest neighbour averaging or polygon contiguity.
}
\value{
An object of class \code{LocalEnv}.
}
\note{
Note that this function is not to interpolate between data points. The calculation of each local environment involves the point itself, so the simple inverse distance function with a power of 2 or more should be used with care.
}
\author{
Seong-Yun Hong
}
\seealso{
\code{\link{SegLocalEnv}}, \code{\link{SegLocalEnv-class}}, \code{\link{spseg}}, \code{dist}
}
\examples{
# Create a random data set with 50 data points and 3 population groups
xy <- matrix(runif(100), ncol = 2)
pop <- matrix(runif(150), ncol = 3)

env1 <- getSegLocalEnv(xy, pop)
summary(env1)

xy.dist <- dist(xy, method = "manhattan")
maxdist <- quantile(as.matrix(xy), 0.75)
env2 <- getSegLocalEnv(xy, pop, maxdist = maxdist, sprel = xy.dist)
summary(env2)

env3 <- getSegLocalEnv(xy, pop, useExp = FALSE, power = 0, maxdist = 0.5)
summary(env3)

env4 <- getSegLocalEnv(xy, pop, useExp = FALSE, maxdist = 0.5)
summary(env4)
}
