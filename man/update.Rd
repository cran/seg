\name{update}
\alias{update}
\alias{update.SegLocalEnv}
\title{Update an Object of Class SegLocalEnv}
\description{
Update an existing object of class \code{SegLocalEnv}.
}
\usage{
  \method{update}{SegLocalEnv}(object, coords, data, env, proj4string, \dots)
}
\arguments{
  \item{object}{an object of class \code{LocalEnv} to be updated.}
  \item{coords}{an optional object of class \code{matrix} containing coordinates of the data points (each row is a point).}
  \item{data}{an optional object of class \code{matrix} containing the population of each group at each data point. The number of rows in \sQuote{data} should equal the number of points in \sQuote{coords}, and the number of columns should be greater than one (at least two population groups are required).}
  \item{env}{an optional object of class \code{matrix} containing the local environment parameters. Must be the same dimensions as \sQuote{data}.}
  \item{proj4string}{an optional projection string of class \code{CRS}.}
  \item{\dots}{ignored.}
}
\value{
An object of class \code{SegLocalEnv}.
}
\author{
Seong-Yun Hong
}
\seealso{
\code{\link{SegLocalEnv}}, \code{\link{getSegLocalEnv}}, \code{\link{SegLocalEnv-class}}
}
\examples{
# Create a random data set with 50 data points and 3 population groups
xy <- matrix(runif(100), ncol = 2)
pop <- matrix(runif(150), ncol = 3)
localenv <- SegLocalEnv(coords = xy, data = pop)
print(localenv)

# Update the projection information
proj <- CRS("+proj=nzmg +datum=nzgd49")
localenv2 <- update(localenv, proj4string = proj)
print(localenv2)
}
