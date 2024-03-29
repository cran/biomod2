% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BIOMOD_ModelingOptions.R
\name{BIOMOD_ModelingOptions}
\alias{BIOMOD_ModelingOptions}
\alias{bm_DefaultModelingOptions}
\title{Configure the modeling options for each selected model}
\usage{
BIOMOD_ModelingOptions(
  GLM = NULL,
  GBM = NULL,
  GAM = NULL,
  CTA = NULL,
  ANN = NULL,
  SRE = NULL,
  FDA = NULL,
  MARS = NULL,
  RF = NULL,
  MAXENT = NULL,
  XGBOOST = NULL
)

bm_DefaultModelingOptions()
}
\arguments{
\item{GLM}{(\emph{optional, default} \code{NULL}) \cr A \code{list} containing GLM options}

\item{GBM}{(\emph{optional, default} \code{NULL}) \cr A \code{list} containing GBM options}

\item{GAM}{(\emph{optional, default} \code{NULL}) \cr A \code{list} containing GAM options}

\item{CTA}{(\emph{optional, default} \code{NULL}) \cr A \code{list} containing CTA options}

\item{ANN}{(\emph{optional, default} \code{NULL}) \cr A \code{list} containing ANN options}

\item{SRE}{(\emph{optional, default} \code{NULL}) \cr A \code{list} containing SRE options}

\item{FDA}{(\emph{optional, default} \code{NULL}) \cr A \code{list} containing FDA options}

\item{MARS}{(\emph{optional, default} \code{NULL}) \cr A \code{list} containing MARS options}

\item{RF}{(\emph{optional, default} \code{NULL}) \cr A \code{list} containing RF options}

\item{MAXENT}{(\emph{optional, default} \code{NULL}) \cr A \code{list} 
containing MAXENT options}

\item{XGBOOST}{(\emph{optional, default} \code{NULL}) \cr A \code{list} 
containing XGBOOST options}
}
\value{
A \code{\link{BIOMOD.models.options}} object that can be used to build species distribution 
model(s) with the \code{\link{BIOMOD_Modeling}} function.
}
\description{
Parametrize and/or tune \pkg{biomod2}'s single models options.
}
\details{
This function allows advanced user to change some default parameters of \pkg{biomod2} inner 
models. \cr 10 single models are available within the package, and their options can be set 
with this function through \code{list} objects.

The \code{\link{bm_DefaultModelingOptions}} function prints all default parameter values for 
all available models. \cr This output can be copied and pasted to be used as is (with wanted 
changes) as function arguments (see \href{https://biomodhub.github.io/biomod2/reference/BIOMOD_ModelingOptions.html#examples}{Examples}).

Below is the detailed list of all modifiable parameters for each available model.
}
\section{GLM }{
 (\code{\link[stats]{glm}})
\itemize{
  \item{\code{myFormula}}{ : a typical \code{formula} object (see 
  \href{https://biomodhub.github.io/biomod2/reference/BIOMOD_ModelingOptions.html#examples}{Examples}). \cr If not \code{NULL}, \code{type} 
  and \code{interaction.level} parameters are switched off. \cr 
  You can choose to either :
  \itemize{
    \item{generate automatically the GLM formula with the following parameters :
    \itemize{
      \item{\code{type = 'quadratic'}}{ : formula given to the model, must 
      be \code{simple}, \code{quadratic} or \code{polynomial}}
      \item{\code{interaction.level = 0}}{ : an \code{integer} corresponding 
      to the interaction level between considered variables considered (\emph{be aware that 
      interactions quickly enlarge the number of effective variables used into the GLM !})}
    }
    }
    \item{or construct specific formula}
  }
  }
  \item{\code{test = 'AIC'}}{ : information criteria for the stepwise 
  selection procedure, must be \code{AIC} (\emph{Akaike Information Criteria}, \code{BIC} 
  (\emph{Bayesian Information Criteria}) or \code{none} (\emph{consider only the full model, 
  no stepwise selection, but this can lead to convergence issue and strange results !})}
  \item{\code{family = binomial(link = 'logit')}}{ : a \code{character} 
  defining the error distribution and link function to be used in the model, mus be a family 
  name, a family function or the result of a call to a family function (see \link{family}) 
  (\emph{so far, \pkg{biomod2} only runs on presence-absence data, so binomial family is the 
  default !})}
  \item{\code{control}}{ : a \code{list} of parameters to control the fitting process (passed to 
  \code{\link{glm.control}})}
}
}

\section{GBM }{
 (default \code{\link[gbm]{gbm}})

\emph{Please refer to \code{\link[gbm]{gbm}} help file for more details.}
\itemize{
  \item{\code{distribution = 'bernoulli'}}
  \item{\code{n.trees = 2500}}
  \item{\code{interaction.depth = 7}}
  \item{\code{n.minobsinnode = 5}}
  \item{\code{shrinkage = 0.001}}
  \item{\code{bag.fraction = 0.5}}
  \item{\code{train.fraction = 1}}
  \item{\code{cv.folds = 3}}
  \item{\code{keep.data = FALSE}}
  \item{\code{verbose = FALSE}}
  \item{\code{perf.method = 'cv'}}
  \item{\code{n.cores = 1}}
}
}

\section{GAM }{
 (\code{\link[gam]{gam}} or \code{\link[mgcv]{gam}})
\itemize{
  \item{\code{algo = 'GAM_gam'}}{ : a \code{character} defining the chosen GAM function, must 
  be \code{GAM_gam} (see \code{\link[gam]{gam}}), \code{GAM_mgcv} (see \code{\link[mgcv]{gam}}) 
  or \code{BAM_mgcv} (see \code{\link[mgcv]{bam}})}
  \item{\code{myFormula}}{ : a typical \code{formula} object (see 
  \href{https://biomodhub.github.io/biomod2/reference/BIOMOD_ModelingOptions.html#examples}{Examples}). \cr If not \code{NULL}, \code{type} 
  and \code{interaction.level} parameters are switched off. \cr 
  You can choose to either :
  \itemize{
    \item{generate automatically the GAM formula with the following parameters :
    \itemize{
      \item{\code{type = 's_smoother'}}{ : the smoother used to generate the formula}
      \item{\code{interaction.level = 0}}{ : an \code{integer} corresponding 
      to the interaction level between considered variables considered (\emph{be aware that 
      interactions quickly enlarge the number of effective variables used into the GLM !})}
    }
    }
    \item{or construct specific formula}
  }
  }
  \item{\code{k = -1}}{a smooth term in a formula argument to gam, must be \code{-1} or 
  \code{4} (see \pkg{gam} \code{\link[gam]{s}} or \pkg{mgcv} \code{\link[mgcv]{s}})}
  \item{\code{family = binomial(link = 'logit')}}{ : a \code{character} defining 
  the error distribution and link function to be used in the model, mus be a family name, a 
  family function or the result of a call to a family function (see \link{family}) 
  (\emph{so far, \pkg{biomod2} only runs on presence-absence data, so binomial family is the 
  default !})}
  \item{\code{control}}{ : a \code{list} of parameters to control the fitting process (passed to 
  \code{\link[mgcv]{gam.control}} or \code{\link[gam]{gam.control}})}
  \item{some options specific to \code{GAM_mgcv} (\emph{ignored if \code{algo = 'GAM_gam'}})}{
  \itemize{
    \item{\code{method = 'GCV.Cp'})}
    \item{\code{optimizer = c('outer','newton')}}
    \item{\code{select = FALSE}}
    \item{\code{knots = NULL}}
    \item{\code{paramPen = NULL}}
  }
  }
}
}

\section{CTA }{
 (\code{\link[rpart]{rpart}})

\emph{Please refer to \code{\link[rpart]{rpart}} help file for more details.}
\itemize{
  \item{\code{method = 'class'}}
  \item{\code{parms = 'default'}}{ : if \code{'default'}, default \pkg{rpart} 
  \code{parms} value are kept}
  \item{\code{cost = NULL}}
  \item{\code{control}}{ : see \code{\link[rpart]{rpart.control}}}
}
}

\section{ANN }{
 (\code{\link[nnet]{nnet}})
\itemize{
  \item{\code{NbCV = 5}}{ : an \code{integer} corresponding to the number of cross-validation 
  repetitions to find best size and decay parameters}
  \item{\code{size = NULL}}{ : an \code{integer} corresponding to the number of units in the 
  hidden layer. If \code{NULL} then size parameter will be optimized by cross-validation based 
  on model AUC (\code{NbCv} cross-validations ; tested size will be the following : 
  \code{c(2, 4, 6, 8)}). It is also possible to give a \code{vector} of size values to be tested, 
  and the one giving the best model AUC will be kept.}
  \item{\code{decay = NULL}}{ : a \code{numeric} corresponding to weight decay. If \code{NULL} 
  then decay parameter will be optimized by cross-validation based on model AUC (\code{NbCv} 
  cross-validations ; tested size will be the following : \code{c(0.001, 0.01, 0.05, 0.1)}). 
  It is also possible to give a \code{vector} of decay values to be tested, and the one giving 
  the best model AUC will be kept.}
  \item{\code{rang = 0.1}}{ : a \code{numeric} corresponding to the initial random weights on 
  \code{[-rang, rang]}}
  \item{\code{maxit = 200}}{ : an \code{integer} corresponding to the maximum number of 
  iterations}
}
}

\section{SRE }{
 (\code{\link{bm_SRE}})
\itemize{
  \item{\code{quant = 0.025}}{ : a \code{numeric} corresponding to the quantile of 
  '\emph{extreme environmental variable}' removed to select species envelops}
}
}

\section{FDA }{
 (\code{\link[mda]{fda}})

\emph{Please refer to \code{\link[mda]{fda}} help file for more details.}
\itemize{
  \item{\code{method = 'mars'}}
  \item{\code{add_args = NULL}}{ : a \code{list} of additional parameters to \code{method} and 
  given to the \code{...} options of \code{\link[mda]{fda}} function}
}
}

\section{MARS }{
 (\code{\link[earth]{earth}})

\emph{Please refer to \code{\link[earth]{earth}} help file for more details.}
\itemize{
  \item{\code{myFormula}}{ : a typical \code{formula} object (see 
  \href{https://biomodhub.github.io/biomod2/reference/BIOMOD_ModelingOptions.html#examples}{Examples}). \cr If not \code{NULL}, \code{type} 
  and \code{interaction.level} parameters are switched off. \cr 
  You can choose to either :
  \itemize{
    \item{generate automatically the MARS formula with the following parameters :
    \itemize{
      \item{\code{type = 'simple'}}{ : formula given to the model, must 
      be \code{simple}, \code{quadratic} or \code{polynomial}}
      \item{\code{interaction.level = 0}}{ : an \code{integer} corresponding 
      to the interaction level between considered variables considered (\emph{be aware that 
      interactions quickly enlarge the number of effective variables used into the MARS !})}
    }
    }
    \item{or construct specific formula}
  }
  }
  \item{\code{nk = NULL}}{ : an \code{integer} corresponding to the maximum number of model 
  terms. \cr 
  If \code{NULL} default MARS function value is used : \code{max(21, 2 * nb_expl_var + 1)}}
  \item{\code{penalty = 2}}
  \item{\code{thresh = 0.001}}
  \item{\code{nprune = NULL}}
  \item{\code{pmethod = 'backward'}}
}
}

\section{RF }{
 (\code{\link[randomForest]{randomForest}})
\itemize{
  \item{\code{do.classif = TRUE}}{ : if \code{TRUE} \emph{random.forest classification} will 
  be computed, otherwise \emph{random.forest regression} will be done}
  \item{\code{ntree = 500}}
  \item{\code{mtry = 'default'}}
  \item{\code{sampsize = NULL}}
  \item{\code{nodesize = 5}}
  \item{\code{maxnodes = NULL}}
}
}

\section{MAXENT }{
 (\url{https://biodiversityinformatics.amnh.org/open_source/maxent/})
\itemize{
  \item{\code{path_to_maxent.jar = getwd()}}{ : a \code{character}
  corresponding to \pkg{maxent.jar} file link} 
  
  \item{\code{memory_allocated = 512}}{ : an \code{integer} corresponding to
  the amount of memory (in Mo) reserved for \code{java} to run
  \code{MAXENT}, must be \code{64}, \code{128}, \code{256},
  \code{512}, \code{1024}... or \code{NULL} to use default \code{java}
  memory limitation parameter}
  
  \item{\code{initial_heap_size = NULL}}{ : a \code{character} initial heap
  space (shared memory space) allocated to java. Argument transmitted to
  \code{-Xms} when calling java. Used in \code{\link{BIOMOD_Projection}} but
  not in \code{\link{BIOMOD_Modeling}}. Values can be \code{1024K},
  \code{4096M}, \code{10G} ... or \code{NULL} to use default \code{java}
  parameter}
  
  \item{\code{max_heap_size = NULL}}{ : a \code{character} initial heap
  space (shared memory space) allocated to java. Argument transmitted to
  \code{-Xmx} when calling java. Used in \code{\link{BIOMOD_Projection}} but
  not in \code{\link{BIOMOD_Modeling}}. Must be larger than
  \code{initial_heap_size}. Values can be \code{1024K}, \code{4096M},
  \code{10G} ... or \code{NULL} to use default \code{java} parameter}
  
  \item{\code{background_data_dir}}{ : a \code{character} corresponding to
  directory path where explanatory variables are stored as \code{ASCII}
  files (raster format). If specified, \code{MAXENT} will generate
  its own background data from explanatory variables rasters (as usually
  done in \code{MAXENT} studies). Otherwise \pkg{biomod2} pseudo-absences
  will be used (see \code{\link{BIOMOD_FormatingData}})}
  
  \item{\code{maximumbackground}}{ : an \code{integer} corresponding to the
  maximum number of background data to sample if the
  \code{background_data_dir} parameter has been set}
  
  \item{\code{maximumiterations = 200}}{ : an \code{integer} corresponding
  to the maximum number of iterations to do} 
  
  \item{\code{visible = FALSE}}{ : a \code{logical} to make the
  \code{MAXENT} user interface available}
  
  \item{\code{linear = TRUE}}{ : a \code{logical} to allow linear features
  to be used} \item{\code{quadratic = TRUE}}{ : a \code{logical} to allow
  quadratic features to be used} 
  
  \item{\code{product = TRUE}}{ : a \code{logical} to allow product features
  to be used}
  
  \item{\code{threshold = TRUE}}{ : a \code{logical} to allow threshold
  features to be used}
  
  \item{\code{hinge = TRUE}}{ : a \code{logical} to allow hinge features to
  be used} \item{\code{lq2lqptthreshold = 80}}{ : an \code{integer}
  corresponding to the number of samples at which product and threshold
  features start being used} 
  
  \item{\code{l2lqthreshold = 10}}{ : an
  \code{integer} corresponding to the number of samples at which quadratic
  features start being used} 
  
  \item{\code{hingethreshold = 15}}{ : an
  \code{integer} corresponding to the number of samples at which hinge
  features start being used}
  
  \item{\code{beta_threshold = -1.0}}{ : a
  \code{numeric} corresponding to the regularization parameter to be applied
  to all threshold features (\emph{negative value enables automatic
  setting})} 
  
  \item{\code{beta_categorical = -1.0}}{ : a \code{numeric}
  corresponding to the regularization parameter to be applied to all
  categorical features (\emph{negative value enables automatic setting})}
  
  \item{\code{beta_lqp = -1.0}}{ : a \code{numeric} corresponding to the
  regularization parameter to be applied to all linear, quadratic and
  product features (\emph{negative value enables automatic setting})}
  
  \item{\code{beta_hinge = -1.0}}{ : a \code{numeric} corresponding to the
  regularization parameter to be applied to all hinge features
  (\emph{negative value enables automatic setting})}
  
  \item{\code{betamultiplier = 1}}{ : a \code{numeric} to multiply all
  automatic regularization parameters \cr (\emph{higher number gives a more
  spread-out distribution})} 
  
  \item{\code{defaultprevalence = 0.5}}{ : a
  \code{numeric} corresponding to the default prevalence of the species \cr
  (\emph{probability of presence at ordinary occurrence points})}
}
}

\section{XGBOOST }{
 (default \code{\link[xgboost]{xgboost}})

\emph{Please refer to \code{\link[xgboost]{xgboost}} help file for more details.}
\itemize{
  \item{\code{max.depth = 5}}
  \item{\code{eta = 0.1}}
  \item{\code{nrounds = 512}}
  \item{\code{objective = "binary:logistic"}}
  \item{\code{nthread = 1}}
  }
}

\examples{
library(terra)

# Load species occurrences (6 species available)
data(DataSpecies)
head(DataSpecies)

# Select the name of the studied species
myRespName <- 'GuloGulo'

# Get corresponding presence/absence data
myResp <- as.numeric(DataSpecies[, myRespName])

# Get corresponding XY coordinates
myRespXY <- DataSpecies[, c('X_WGS84', 'Y_WGS84')]

# Load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
data(bioclim_current)
myExpl <- terra::rast(bioclim_current)

\dontshow{
myExtent <- terra::ext(0,30,45,70)
myExpl <- terra::crop(myExpl, myExtent)
}

# ---------------------------------------------------------------#
# Print default modeling options
bm_DefaultModelingOptions()

# Create default modeling options
myBiomodOptions <- BIOMOD_ModelingOptions()
myBiomodOptions

# # Part (or totality) of the print can be copied and customized
# # Below is an example to compute quadratic GLM and select best model with 'BIC' criterium
# myBiomodOptions <- BIOMOD_ModelingOptions(
#   GLM = list(type = 'quadratic',
#              interaction.level = 0,
#              myFormula = NULL,
#              test = 'BIC',
#              family = 'binomial',
#              control = glm.control(epsilon = 1e-08,
#                                    maxit = 1000,
#                                    trace = FALSE)))
# myBiomodOptions
# 
# # It is also possible to give a specific GLM formula
# myForm <- 'Sp277 ~ bio3 + log(bio10) + poly(bio16, 2) + bio19 + bio3:bio19'
# myBiomodOptions <- BIOMOD_ModelingOptions(GLM = list(myFormula = formula(myForm)))
# myBiomodOptions


}
\seealso{
\code{\link{BIOMOD_Tuning}}, \code{\link{BIOMOD_Modeling}}

Other Main functions: 
\code{\link{BIOMOD_EnsembleForecasting}()},
\code{\link{BIOMOD_EnsembleModeling}()},
\code{\link{BIOMOD_FormatingData}()},
\code{\link{BIOMOD_LoadModels}()},
\code{\link{BIOMOD_Modeling}()},
\code{\link{BIOMOD_PresenceOnly}()},
\code{\link{BIOMOD_Projection}()},
\code{\link{BIOMOD_RangeSize}()},
\code{\link{BIOMOD_Tuning}()}
}
\author{
Damien Georges, Wilfried Thuiller
}
\concept{Main functions}
\keyword{models}
\keyword{options}
