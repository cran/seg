\name{SegSpatial-class}
\docType{class}
\alias{SegSpatial-class}
\alias{SegSpatialExt-class}
\alias{show,SegSpatial-method}
\alias{show,SegSpatialExt-method}
\alias{coerce,SegSpatial,list-method}
\alias{coerce,list,SegSpatial-method}
\alias{coerce,SegSpatialExt,list-method}
\alias{print.SegSpatial-method}
\alias{as.list.SegSpatial-method}
\alias{plot.SegSpatialExt-method}
\alias{print.SegSpatialExt-method}
\alias{as.list.SegSpatialExt-method}

\title{Class SegSpatial}
\description{Class for spatial segregation index values.}
\section{Objects from the Class}{
Objects can be created by the function \code{\link{spseg}}, or \code{\link{SegSpatial}}.
}
\section{Slots}{
  \describe{
    \item{d}{an object of class \code{numeric}, indicating the degree of the spatial dissimilarity index.}
    \item{r}{an object of class \code{numeric}, indicating the degree of the spatial diversity index.}
    \item{h}{an object of class \code{numeric}, indicating the degree of the spatial information theory index.}
    \item{p}{an object of class \code{matrix}, showing the spatial exposure/isolation of all population groups.}
  }
}
\section{Methods}{
  \describe{
    \item{coerce}{\code{signature(from = "SegSpatial", to = "list")}: coerce an object of class \code{SegSpatial} to an object of class \code{list}.}
    \item{coerce}{\code{signature(from = "list", to = "SegSpatial")}: coerce an object of class \code{list} to an object of class \code{SegSpatial}.}
    \item{as.list}{\code{signature(x = "SegSpatial")}: coerce an object of class \code{SegSpatial} to an object of class \code{list}.}
    \item{show}{\code{signature(object = "SegSpatial")}: print the spatial segregation values.}
    \item{print}{\code{signature(x = "SegSpatial")}: same as \code{show}.}
  }
}
\author{Seong-Yun Hong}

\seealso{
\code{\link{spseg}}
}
\keyword{classes}
