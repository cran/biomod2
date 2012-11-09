\name{BIOMOD.EnsembleModeling.out-method}

\alias{getEMalgos}
\alias{getEMalgos,BIOMOD.EnsembleModeling.out-method}
\alias{getEMbuiltModels}
\alias{getEMbuiltModels,BIOMOD.EnsembleModeling.out-method}
\alias{getEMeval}
\alias{getEMeval,BIOMOD.EnsembleModeling.out-method}
\alias{getEMkeptModels}
\alias{getEMkeptModels,BIOMOD.EnsembleModeling.out-method}

\title{ BIOMOD.EnsembleModeling.out getters }

\description{ 
Functions to get attributs of \code{\link[biomod2]{BIOMOD_EnsembleModeling}} outputs
}

\usage{
  getEMalgos(obj, model)
  getEMbuiltModels(obj)
  getEMeval(obj, ...)
  getEMkeptModels(obj, model)
}



\arguments{
  \item{obj}{ \code{"\link[=BIOMOD.EnsembleModeling.out-class]{BIOMOD.EnsembleModeling.out}"} object }
  \item{model}{"character" or "numeric" refering to model names (getEMbuiltModels()) or model id }
  \item{...}{ extra arguments (see details) }
}

\details{

  \ldots{ } available values : 
  
  \itemize{
  
    \item{\bold{getEMeval}}{
      \itemize{
        \item{\code{model}:}{"character", model name }
        \item{\code{met}:}{"character", evaluation metric name}
      }
    }
  } 
}


\value{

  \enumerate{
    \item{\bold{getEMalgos: }}{ a \code{character} vector indicating set of ensemble-modeling algorithms ran}
    
    \item{\bold{getEMbuiltModels: }}{ a \code{character} vector indicating the names of ensemble models computed }
    
    \item{\bold{getEMeval: }}{ an \code{array} or a \code{data.frame} containing ensemble models evaluation scores}
    
    \item{\bold{getEMkeptModels: }}{ an \code{character} vector indicating names of selected models for ensemble-models building}
    
  }
}

\author{ Wilfried Thuiller, Damien Georges}

\seealso{ \code{\link{help}}}

\keyword{ models }