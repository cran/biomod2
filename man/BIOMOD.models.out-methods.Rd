\name{BIOMOD.models.out-methods}

\alias{getModelsPrediction}
\alias{getModelsPrediction,BIOMOD.models.out-method}
\alias{getModelsPredictionEval}
\alias{getModelsPredictionEval,BIOMOD.models.out-method}
\alias{getModelsEvaluations}
\alias{getModelsEvaluations,BIOMOD.models.out-method}
\alias{getModelsVarImport}
\alias{getModelsVarImport,BIOMOD.models.out-method}
\alias{getModelsOptions}
\alias{getModelsOptions,BIOMOD.models.out-method}
\alias{getModelsInputData}
\alias{getModelsInputData,BIOMOD.models.out-method}
\alias{getModelsBuiltModels}
\alias{getModelsBuiltModels,BIOMOD.models.out-method}

\title{ BIOMOD.models.out getters }

\description{ 
Functions to get attributs of \code{\link[biomod2]{BIOMOD_Modeling}} outputs
}

\usage{
  getModelsPrediction(obj, ...)
  getModelsPredictionEval(obj, ...)
  getModelsEvaluations(obj)
  getModelsVarImport(obj)
  getModelsOptions(obj)
  getModelsInputData(obj, ...)
  getModelsBuiltModels(obj)
}



\arguments{
  \item{obj}{ \code{"\link[=BIOMOD.models.out-class]{BIOMOD.models.out}"} object }
  \item{...}{ extra arguments (see details) }
}

\details{

  \ldots{ } available values : 
  
  \itemize{
  
    \item{\bold{getModelsPrediction & getModelsPredictionEval}}{
      \itemize{
        \item{\code{as.data.frame}:}{logical. If TRUE, models predictions will be returned as \code{data.frame} rather than \code{array}}
      }
    }
    
    \item{\bold{getModelsInputData}}{
      \itemize{
        \item{\code{subinfo}:}{character. Flag defining a specific information to extract from \code{"\link[=BIOMOD.formated.data-class]{BIOMOD.formated.data}"} object. Supported values are: }
          \itemize{
            \item{\code{NULL}:}{ (default) the whole \code{"\link[=BIOMOD.formated.data-class]{BIOMOD.formated.data}"} object is returned }
            \item{\code{'MinMax'}:}{ All explanatory variables ranges returned}
            \item{\code{'expl.var'}:}{ Explanatory variables \code{data.frame} returned}
          }
      }
    }
  }
}


\value{

  \enumerate{
    \item{\bold{getModelsPrediction: }}{ an \code{array} (or a \code{data.frame}) containing models predictions over calibrating and testing data (those used for evaluate models)}
    
    \item{\bold{getModelsPredictionEval: }}{ an \code{array} (or a \code{data.frame}) containing models predictions over evaluating data (if referenced / those used for evaluate models) }
    
    \item{\bold{getModelsEvaluations: }}{ an \code{array} containing models variables importances}
    
    \item{\bold{getModelsVarImport: }}{ an \code{array} containing models variables importances}
    
    \item{\bold{getModelsOptions: }}{ a \code{"\link[=BIOMOD.Model.Options-class]{BIOMOD.Model.Options}"} reporting options used to build individual models}
    
    \item{\bold{getModelsInputData: }}{a \code{"\link[=BIOMOD.formated.data-class]{BIOMOD.formated.data}"} object containing data used for models building and evavaluation, or a part of this object}
    
    \item{\bold{getModelsBuiltModels: }}{a character vector giving the names of models succefully computed} 
  
  }
}

\author{ Wilfried Thuiller, Damien Georges}

\seealso{ \code{\link{help}}}

\keyword{ models }