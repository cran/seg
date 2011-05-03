\name{SegLocalEnv-class}
\docType{class}
\alias{SegLocalEnv-class}
\alias{show,SegLocalEnv-method}
\alias{coerce,SegLocalEnv,list-method}
\alias{coerce,list,SegLocalEnv-method}
\alias{coerce,SegLocalEnv,SpatialPointsDataFrame-method}
\alias{coerce,SpatialPointsDataFrame,SegLocalEnv-method}
\alias{coerce,SpatialPolygonsDataFrame,SegLocalEnv-method}
\alias{print.SegLocalEnv-method}
\alias{summary.SegLocalEnv-method}
\alias{as.list.SegLocalEnv-method}
\alias{plot.SegLocalEnv-method}
\alias{points.SegLocalEnv-method}

\title{Class SegLocalEnv}
\description{Class for local environments of spatial data points.}
\section{Objects from the Class}{
Objects can be created by the function \code{SegLocalEnv}.
}
\section{Slots}{
  \describe{
    \item{coords}{an object of class \code{matrix}, containing the coordinates (each row is a point).}
    \item{data}{an object of class \code{matrix}, containing the population composition at each point (each column is a population group).}
    \item{env}{an object of class \code{matrix}, containing the population composition of the local environment for each point.}
    \item{proj4string}{an object of class \code{CRS}, projection string.}
  }
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "SegLocalEnv", to = "list")}: coerce an object of class \code{SegLocalEnv} to an object of class \code{list}.}
    \item{coerce}{\code{signature(from = "SegLocalEnv", to = "SpatialPointsDataFrame")}: coerce an object of class \code{SegLocalEnv} to an object of class \code{SpatialPointsDataFrame}.}
    \item{coerce}{\code{signature(from = "list", to = "SegLocalEnv")}: coerce an object of class \code{list} to an object of class \code{SegLocalEnv}.}
    \item{coerce}{\code{signature(from = "SpatialPointsDataFrame", to = "SegLocalEnv")}: coerce an object of class \code{SpatialPointsDataFrame} to an object of class \code{SegLocalEnv}.}
    \item{coerce}{\code{signature(from = "SpatialPolygonsDataFrame", to = "SegLocalEnv")}: coerce an object of class \code{SpatialPolygonsDataFrame} to an object of class \code{SegLocalEnv}.}
    \item{as.list}{\code{signature(x = "SegLocalEnv")}: coerce an object of class \code{SegLocalEnv} to an object of class \code{list}.}
    \item{show}{\code{signature(object = "SegLocalEnv")}: show the number of points and data columns in an object of class \code{SegLocalEnv}.}
    \item{print}{\code{signature(x = "SegLocalEnv")}: same as \code{show}.}
    \item{summary}{\code{signature(object = "SegLocalEnv")}: summarise the population compositions of points and local environments in an object of class \code{SegLocalEnv}.}
    \item{plot}{\code{signature(x = "SegLocalEnv")}: draw a plot, or plots, of points in an object of class \code{SegLocalEnv}. Use an optional argument \sQuote{which.col} to specify a column of the data that determines the points' sizes. See \sQuote{Examples} below for demonstration.}
    \item{points}{\code{signature(x = "SegLocalEnv")}: draw points in an object of class \code{SegLocalEnv} on an active graphic device.}
  }
}
\author{Seong-Yun Hong}
\seealso{
\code{\link{SegLocalEnv}}, \code{\link{getSegLocalEnv}}, \code{\link{update}}
}
\examples{
# Create a random data set with 50 data points and 3 population groups
xy <- matrix(runif(100), ncol = 2)
pop <- matrix(runif(150), ncol = 3)
localenv <- SegLocalEnv(coords = xy, data = pop)

# Generic functions
print(localenv)
summary(localenv)
par(mfrow = c(1, 3))
plot(localenv, xlab = "x", ylab = "y")
par(mfrow = c(1, 2))
plot(localenv, xlab = "x", ylab = "y", which.col = c(1, 3))
as.list(localenv)
}
\keyword{classes}
