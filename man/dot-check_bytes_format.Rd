% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biomod2_internal.R
\name{.check_bytes_format}
\alias{.check_bytes_format}
\title{Check bytes formatting}
\usage{
.check_bytes_format(this_test, x, varname)
}
\arguments{
\item{x}{string to be transformed}
}
\value{
a boolean
}
\description{
Internal function that check a character string to match a byte
format for Java. e.g. 1024M, 1024m, 1024k or 1024K
}
\keyword{internal}
