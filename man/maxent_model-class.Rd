\name{maxent_model-class}

\docType{class}

\alias{maxent_model-class}
\alias{maxent_model}

% BIOMOD.models.out-method
\alias{show,maxent_model-method}
\alias{predict,maxent_model-method}
\alias{predict.maxent_model.data.frame}
\alias{predict.maxent_model.RasterStack}

\title{ biomod2 MAXENT models objects class and functions}

\description{
  This object represent a MAXENT models as it was a 'classical' \R model. It is produce during \code{\link[biomod2]{BIOMOD_Modeling}} step. It can be printed. \code{\link{predict}} function is also implemented for this object and supports \code{\link{data.frame}} and \code{{RasterStack}} inputs.

}


\section{Slots}{

  \describe{
    \item{\code{sp_name}:}{ "character", species name }
    \item{\code{model_name}:}{ "character", model name }
    \item{\code{model_output_dir}:}{"character", path to model output directory}
    \item{\code{model_options}:}{ "list", list of options used for model building }
  }


}


\author{ Damien Georges }

\examples{
showClass("maxent_model")
}

\keyword{models, maxent, predict}
