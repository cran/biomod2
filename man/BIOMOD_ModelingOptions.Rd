\name{BIOMOD_ModelingOptions}
\alias{BIOMOD_ModelingOptions}

\title{ Configure the modeling options for each selected model}
\description{
  Function to set the different options for each modeling technique.
}


\usage{
BIOMOD_ModelingOptions( GLM = NULL,
                        GBM = NULL,
                        GAM = NULL,
                        CTA = NULL,
                        ANN = NULL,
                        SRE = NULL,
                        FDA = NULL,
                        MARS = NULL,
                        RF = NULL,
                        MAXENT = NULL)
}

\arguments{
\item{GLM}{list of GLM options}
\item{GBM}{list of GBM options}
\item{GAM}{list of GAM options}
\item{CTA}{list of CTA options}
\item{ANN}{list of ANN options}
\item{SRE}{list of SRE options}
\item{FDA}{list of FDA options}
\item{MARS}{list of MARS options}
\item{RF}{list of RF options}
\item{MAXENT}{list of MAXENT options}
}

\details{
The aim of this function is to allow advanced user to change some default parameters of BIOMOD inner models.  
For each modeling technique, options can be set up.  

Each argument have to be put in a list object.

The best way to use this function is to print defaut models options (\code{\link{Print_Default_ModelingOptions}}) or create a default 'BIOMOD.model.option object' and print it in your console. Then copy the output, change only the required parameters, and paste it as function arguments. (see example)  

Here the detailed list of modifiable parameters. They correspond to the traditional parameters that could be setted out for each modeling technique (e.g. ?GLM) 
}

\section{ =-=-= GLM =-=-= (\code{\link[stats]{glm}})}{

  \itemize{
  
    \item{\code{myFormula} : a typical formula object (see example). If not NULL, type and interaction.level args are switched off. 
    You can choose to either: 
    \itemize{
    \item{generate automatically the GLM formula by using the type and interaction.level arguments 
    type : formula given to the model ('simple', 'quadratic' or 'polynomial'). 
    interaction.level : integer corresponding to the interaction level between variables considered. Consider that interactions quickly enlarge the number of effective variables used into the GLM.}
    \item{or construct specific formula}
    }}
    
    \item{\code{test} : Information criteria for the stepwise selection procedure: AIC for Akaike Information Criteria, and BIC for Bayesian Information Criteria ('AIC' or 'BIC'). 'none' is also a supported value which implies to concider only the full model (no stepwise selection). This can lead to convergence issu and strange results.}
    
    \item{\code{family} : a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \link{family} for details of family functions.) . BIOMOD only runs on presence-absence data so far, so binomial family by default.}
    
    \item{\code{control} : a list of parameters for controlling the fitting process. For glm.fit this is passed to \code{\link{glm.control}}.}

  }
}

\section{=-=-= GBM =-=-= (\code{\link[gbm]{gbm}})}{

Please refer to \code{\link[gbm]{gbm}} help file to get the meaning of this options.
  \itemize{
    \item{ n.trees}
    \item{ cv.folds}
    \item{ distribution}
    \item{ interaction.depth}
    \item{ shrinkage}
    \item{ bag.fraction}
    \item{ train.fraction}
  }

}

\section{=-=-= GAM =-=-=  (\code{\link[gam]{gam}} or \code{\link[mgcv]{gam}})}{
  \itemize{
  
    \item{algo : either "GAM_mgcv" (default), "GAM_gam" or "BAM_mgcv" defining the chosen GAM function (see \code{\link[mgcv]{gam}}, \code{\link[gam]{gam}} resp. \code{\link[mgcv]{bam}} for more details)}
    
    \item{\code{myFormula} : a typical formula object (see example). If not NULL, type and interaction.level args are switched off. 
    You can choose to either: 
    \itemize{
    \item{generate automatically the GLM formula by using the type and interaction.level arguments 
    type : the smother used to generate the formula. Only "s_smoother" available at time. 
    interaction.level : integer corresponding to the interaction level between variables considered. Consider that interactions quickly enlarge the number of effective variables used into the GAM. Interaction are not considered if you choosed "GAM_gam" algo}
    \item{or construct specific formula}
    }}
    
    \item{k : a smooth term in a formula argument to gam (see \pkg{gam} \code{\link[gam]{s}} or \pkg{mgcv} \code{\link[mgcv]{s}})}
%    \item{test : Information criteria for the stepwise selection procedure: AIC for Akaike Information Criteria, and BIC for Bayesian Information Criteria ('AIC' or 'BIC')}  
    \item{family : a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \link{family} for details of family functions.) . BIOMOD only runs on presence-absence data so far, so binomial family by default. }
    \item{control : see \code{\link[mgcv]{gam.control}} or \code{\link[gam]{gam.control}}}  
  }

}

\section{=-=-= CTA =-=-=  (\code{\link[rpart]{rpart}})}{

  Please refer to \code{\link[rpart]{rpart}} help file to get the meaning of the following options.
  \itemize{
    \item{method ('class')}
    \item{parms ('default')}
    \item{cost}
    \item{control}
  }

NOTE: for method and parms, you can give a 'real' value as described in the rpart help file or 'default' that implies default \code{\link[rpart]{rpart}} values.
}

\section{=-=-= ANN =-=-=  (\code{\link[nnet]{nnet}})}{

  \itemize{
    \item{NbCV : nb of cross validation to find best size and decay parameters}
    \item{rang : Initial random weights on [-rang, rang]}
    \item{maxit : maximum number of iterations. Default 100}
  }

}

\section{=-=-= SRE =-=-= }{
  \itemize{
    \item{quant : quantile of 'extreme environmental variable' removed for selection of species envelops}
  }

}

\section{=-=-= FDA =-=-=  (\code{\link[mda]{fda}})}{

  Please refer to \code{\link[mda]{fda}} help file to get the meaning of these options.
  \itemize{
    \item{method}
  }
}

\section{=-=-= MARS -=-= (\code{\link[mda]{mars}})}{

  Please refer to \code{\link[mda]{mars}} help file to get the meaning of these options.
  \itemize{
    \item{degree}
    \item{penalty}
    \item{thresh}
    \item{prune}
  }

}

\section{=-=-= RF -=-=-= (\code{\link[randomForest]{randomForest}})}{

  \itemize{
    \item{do.classif : if TRUE classification random.forest computed else regression random.forest will be done}
    \item{ntree}
    \item{mtry ('default)}
  }

  NOTE: for mtry, you can give a 'real' value as described in randomForest help file or 'default' that implies default randomForest values
  
}

\section{=-=-= MAXENT -=  \url{http://www.cs.princeton.edu/~schapire/maxent/}}{
  \itemize{
    \item{path_to_maxent.jar : character, the link to \pkg{maxent.jar} file (the working directory by default) }
    \item{maximumiterations : integer (default 200), maximum iteration done}
    \item{visible : logical (default FALSE), make the Maxent user interface visible}
    \item{linear : logical (default TRUE), allow linear features to be used}
    \item{quadratic : logical (default TRUE), allow quadratic features to be used}
    \item{product : logical (default TRUE), allow product features to be used}
    \item{threshold : logical (default TRUE), allow threshold features to be used}
    \item{hinge : logical (default TRUE), allow hinge features to be used}
    \item{lq2lqptthreshold : integer (default 80), number of samples at which product and threshold features start being used}
    \item{l2lqthreshold : integer (default 10), number of samples at which quadratic features start being used}
    \item{hingethreshold : integer (default 15), number of samples at which hinge features start being used}
    \item{beta_threshold : numeric (default -1.0), regularization parameter to be applied to all threshold features; negative value enables automatic setting}
    \item{beta_categorical : numeric (default -1.0), regularization parameter to be applied to all categorical features; negative value enables automatic setting}
    \item{beta_lqp : numeric (default -1.0), regularization parameter to be applied to all linear, quadratic and product features; negative value enables automatic setting}
    \item{beta_hinge : numeric (default -1.0), regularization parameter to be applied to all hinge features; negative value enables automatic setting}
    \item{defaultprevalence : numeric (default 0.5), default prevalence of the species: probability of presence at ordinary occurrence points}
  }

}


\value{ 
 A \code{"\link[=BIOMOD.Model.Options-class]{BIOMOD.Model.Options}"} object given to \code{\link[biomod2]{BIOMOD_Modeling}}
}

\author{ Wilfried Thuiller, Damien Georges }

\examples{
# default BIOMOD.model.option object
myBiomodOptions <- BIOMOD_ModelingOptions()

# print the object
myBiomodOptions

# you can copy a part of the print, change it and custom your options 
# here we want to compute quadratic GLM and select best model with 'BIC' criterium
myBiomodOptions <- BIOMOD_ModelingOptions(
                      GLM = list( type = 'quadratic',
                                  interaction.level = 0,
                                  myFormula = NULL,
                                  test = 'BIC',
                                  family = 'binomial',
                                  control = glm.control(epsilon = 1e-08, maxit = 1000, trace = FALSE) ))
                                  
# check changes was done
myBiomodOptions

# you can prefer to establish your own GLM formula
myBiomodOptions <- BIOMOD_ModelingOptions(
                        GLM = list(myFormula = formula("Sp277 ~ bio3 + log(bio10) + poly(bio16,2) + bio19 + bio3:bio19")))

# check changes was done
myBiomodOptions

# you also can directly print default parameters and then follow the same processus
Print_Default_ModelingOptions()

}

\keyword{ models }
\keyword{ options }