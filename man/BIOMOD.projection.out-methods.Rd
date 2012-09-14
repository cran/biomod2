\name{BIOMOD.projection.out-method}

\alias{getProjection}
\alias{getProjection,BIOMOD.projection.out-method}
\alias{free}
\alias{free,BIOMOD.projection.out-method}

\title{ BIOMOD.projection.out getters }

\description{ 
Functions to get attributs of \code{\link[biomod2]{BIOMOD_Projection}} outputs
}

\usage{
  getProjection(obj, ...)
  free(obj)
}



\arguments{
  \item{obj}{ \code{"\link[=BIOMOD.projection.out-class]{BIOMOD.projection.out}"} object }
  \item{...}{ extra arguments (see details) }
}

\details{

  \ldots{ } available values : 
  
  \itemize{
    \item{\code{model}:}{NULL or "character", defining models you want to extract projections }
    \item{\code{as.data.frame}:}{"logical", return projections into a \code{data.frame} (TRUE) or an \code{array} (FALSE, default)}
  }
}


\value{

  \enumerate{
    \item{\bold{getProjection: }}{ an \code{array} or a \code{data.frame} containing models projections}
    \item{\bold{free: }}{remove projection value from R memory and only keep the path to correspunding file on hard drive.}
  }
}

\author{ Wilfried Thuiller, Damien Georges}

\seealso{ \code{\link{help}}}

\keyword{ models, projection, project, forecast }