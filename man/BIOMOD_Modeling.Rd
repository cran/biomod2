\name{BIOMOD_Modeling}
\alias{BIOMOD_Modeling}

\title{ Run a range of species distribution models}
\description{
  This function allows to calibrate and evaluate a range of species distribution models techniques run over a given species. Calibrations are made on the whole sample or a random subpart. The predictive power of the different models is estimated using a range of evaluation metrics.
}
\usage{
BIOMOD_Modeling( data, 
                 models = c('GLM','GBM','GAM','CTA','ANN',
                            'SRE','FDA','MARS','RF','MAXENT'), 
                 models.options = NULL, 
                 NbRunEval=1, 
                 DataSplit=100, 
                 Yweights=NULL,
                 Prevalence=NULL,
                 VarImport=0, 
                 models.eval.meth = c('KAPPA','TSS','ROC'), 
                 SaveObj = TRUE,
                 rescal.all.models = TRUE,
                 do.full.models = TRUE)
}

\arguments{
\item{data}{\code{BIOMOD.formated.data} object returned by \code{\link[biomod2]{BIOMOD_FormatingData}}}
\item{models}{vector of models names choosen among 'GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF' and 'MAXENT' } 
\item{models.options}{\code{BIOMOD.models.options} object returned by \code{\link[biomod2]{BIOMOD_ModelingOptions}}}
\item{NbRunEval}{Number of Evaluation run} 
\item{DataSplit}{\% of data used to calibrate the models, the remaining part will be used for testing }
\item{Yweights}{response points weights}
\item{Prevalence}{either \code{NULL} (default) or a 0-1 numeric used to build 'weighted response weights'}
\item{VarImport}{Number of permutation to estimate variable importance}
\item{models.eval.meth}{vector of names of evaluation metric among 'KAPPA', 'TSS', 'ROC', 'FAR', 'SR', 'ACCURACY', 'BIAS', 'POD', 'CSI' and 'ETS'}
\item{SaveObj}{keep all results and outputs on hard drive or not (NOTE: strongly recommended)}
\item{rescal.all.models}{if true, all model prediction will be rescaled with a binomial GLM}
\item{do.full.models}{if true, models calibrated and evaluated with the whole dataset are done }
}

\details{

\enumerate{

  \item{\bold{data}}{
  
  If you have decide to add pseudo absences to your original dataset (see \code{\link[biomod2]{BIOMOD_FormatingData}}), NbPseudoAbsences * \code{NbRunEval + 1} models will be created. 
  
  }
  
  \item{\bold{models}}{
  
    The set of models to be calibrated on the data. 10 modeling techniques are currently available: 
    \itemize{
      \item{GLM : Generalized Linear Model (\code{\link[stats]{glm}})}
      \item{GAM : Generalized Additive Model (\code{\link[gam]{gam}}, \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}}, see \code{\link[biomod2]{BIOMOD_ModelingOptions} for details on algorithm selection})}
      \item{GBM : Generalized Boosting Model or usually called Boosted Regression Trees (\code{\link[gbm]{gbm}})}
      \item{CTA : Classification Tree Analysis (\code{\link[rpart]{rpart}})}
      \item{ANN : Artificial Neural Network (\code{\link[nnet]{nnet}})}
      \item{SRE : Surface Range Envelop or usually called BIOCLIM}
      \item{FDA : Flexible Discriminant Analysis (\code{\link[mda]{fda}})}
      \item{MARS : Multiple Adaptive Regression Splines (\code{\link[mda]{mars}})}
      \item{RF : Random Forest (\code{\link[randomForest]{randomForest}})}
      \item{MAXENT : Maximum Entropy (\url{http://www.cs.princeton.edu/~schapire/maxent/})}
    }
    
  }
  

   \item{\bold{NbRunEval & DataSplit}}{
   
     As already explained in the \code{\link{BIOMOD_FormatingData}} help file, the common trend is to split the original dataset into two subsets, one to calibrate the models, and another one to evaluate them. Here we provide the possibility to repeat this process (calibration and evaluation) N times (\code{NbRunEval} times). The proportion of data kept for calibration is determined by the \code{DataSplit} argument (100\% - \code{DataSplit} will be used to evaluate the model). This sort of cross-validation allows to have a quite robust test of the models when independent data are not available. Each technique will also be calibrated on the complete original data. All the models produced by BIOMOD and their related informations are saved on the hard drive. 
     
   }
  
  \item{\bold{Yweights & Prevalence}}{
  
    Allows to give more or less weight to some particular observations. If this argument is kept to NULL, each site (presence or absence) has the same weight. In the particular case that pseudo-absence data have been generated \code{\link{BIOMOD_FormatingData}}, weights are calculated in a way to obtain the given prevalence (0.5 if \code{NULL} meaning that the presences will have the same importance than the absences in the calibration process of the models). Automatically created \code{Yweights} will be composed by integers to prevent from some modelling issues.
    Note that \code{Prevalence} argument will always be ignored if \code{Yweights} are defined.
    
    
  }
  
  \item{\bold{models.eval.meth}}{
  
    The available evaluations methods are : 
    \itemize{
      \item{\sQuote{ROC} : Relative Operating Characteristic}
      \item{\sQuote{KAPPA} : Cohen's Kappa (Heidke skill score) }
      \item{\sQuote{TSS} : True kill statistic (Hanssen and Kuipers discriminant, Peirce's skill score)}
      \item{\sQuote{FAR} : False alarm ratio }
      \item{\sQuote{SR} : Success ratio}
      \item{\sQuote{ACCURANCY} : Accuracy (fraction correct)}
      \item{\sQuote{BIAS} : Bias score (frequency bias)}
      \item{\sQuote{POD} : Probability of detection (hit rate)}
      \item{\sQuote{CSI} : Critical success index (threat score)}
      \item{\sQuote{ETS} : Equitable threat score (Gilbert skill score)}
    }
    Some of them are rescaled to have all an optimum at 1. You can choose one of more (vector) evaluation metric. By Default, only 'KAPPA', 'TSS' and 'ROC' evaluation are done. Please refer to the CAWRC website (\url{http://www.cawcr.gov.au/projects/verification/#Methods_for_dichotomous_forecasts}) to get detailled description of each metric.
    
  }
  
  \item{\bold{SaveObj}}{
  
  If this argument is set to False, it may prevent the evaluation of the \sQuote{ensemble modelled} models in further steps. We strongly recommend to always keep this argument \code{TRUE} even it asks for free space onto the hard drive.
  }
  
  \item{\bold{rescal.all.models}}{
  
  Some categorical models have to be rescaled in every case (\sQuote{MARS}, \sQuote{FDA}, \sQuote{ANN}). But It may be interesting to rescale all model computed to ensure that they will produced comparable predictions (0-1000 ladder). That's particularly useful when you do some ensemble forecasting to remove the scale prediction effect (the more extended projections are, the more they influence ensemble forecasting results). 
  }

  \item{\bold{do.full.models}}{
  
  Building models with all information available may be usefull in some particular cases (i.e. rare species with few presences points). The main drawback of this method is that, if you don't give separated data for models evaluation, your models will be evaluated with the same data that the ones used for calibration. Thats will lead to over-optimistic evaluation scores. Be carefull whith this '_Full' models interpretation.
  }
} % end of enumerate
} % end of detail 

\value{
  A BIOMOD.models.out object
    
  See \code{"\link[=BIOMOD.models.out-class]{BIOMOD.models.out}"} for details.
  
  Additional objects are stored out of R in two different directories for memory storage purposes. They are created by the function directly on the root
  of your working directory set in R ("models" directory). This one contains each calibrated model for each repetition and pseudo-absence run. A hidden folder \file{.DATA_BIOMOD} contains some files (predictions, original dataset copy, pseudo absences chosen...) used by other functions like \code{\link[biomod2]{BIOMOD_Projection}} or \code{\link[biomod2]{BIOMOD_EnsembleModeling}} .
  
  The models are currently stored as objects to be read exclusively in R. To load them back (the same stands for all objects stored on the hard disk) 
  use the \code{\link{load}} function (see examples section below). 
   
}


\author{ 
Wilfried Thuiller, Damien Georges, Robin Engler
}

\seealso{
\code{\link{BIOMOD_FormatingData}},  \code{\link{BIOMOD_ModelingOptions}}, \code{\link{BIOMOD_Projection}}
}

\examples{

# 0. Load data & Selecting Data

# species occurances
species_occ <- read.csv(system.file("external/species/species_occ.csv",package="biomod2"))

# we consider only presences of MyocastorCoypus species
myRespName <- 'MyocastorCoypus'
myRespCoord <- species_occ[which(!is.na(species_occ[,myRespName])),c('x','y')]
myResp <- as.numeric(na.omit(species_occ[,myRespName]))

# Environemental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
myExpl = raster::stack(system.file("external/climat/current/bio3.grd",package="biomod2"),
                       system.file("external/climat/current/bio4.grd",package="biomod2"),
                       system.file("external/climat/current/bio7.grd",package="biomod2"),
                       system.file("external/climat/current/bio11.grd",package="biomod2"),
                       system.file("external/climat/current/bio12.grd",package="biomod2"))

# 1. Formating Data
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespCoord,
                                     resp.name = myRespName,
                                     PA.nb.rep = 2,
                                     PA.nb.absences = 200,
                                     PA.strategy = 'random')
                                                                     
# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

# 3. Doing Modelisation

myBiomomodModelOut <- BIOMOD_Modeling( myBiomodData, 
                                       models = c('SRE','CTA','RF'), 
                                       models.options = myBiomodOption, 
                                       NbRunEval=1, 
                                       DataSplit=80, 
                                       Yweights=NULL, 
                                       VarImport=3, 
                                       models.eval.meth = c('TSS','ROC'),
                                       SaveObj = TRUE )
                                       
## print a summary of modeling stuff
myBiomomodModelOut
                                       
                                    
}

\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ tree }
