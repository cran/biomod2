\name{Find.Optim.Stat}
\alias{Find.Optim.Stat}

\title{ Find the threshold driving to the best evaluation score }

\description{ \code{Find.Optim.Stat} is an internal \pkg{biomod2} function to find the threshold to convert continuous values into binary ones leading to the best score for a given evaluation metric.
}

\usage{
Find.Optim.Stat(Stat='TSS',
                Fit,
                Obs,
                Precision = 5,
                Fixed.thresh = NULL)
}

\arguments{
  \item{Stat}{  either 'ROC', TSS', 'KAPPA', 'ACCURACY', 'BIAS', 'POD', 'FAR', 'POFD', 'SR', 'CSI', 'ETS', 'HK', 'HSS', 'OR' or 'ORSS'}
  \item{Fit}{ vector of fitted values (continuous)}
  \item{Obs}{ vector of observed values (binary)}
  \item{Precision}{ threshold test steps (a threshold will be test each \code{Precision})}
  \item{Fixed.thresh}{ if not \code{NULL}, the only threshold value tested }
}

\details{
  Please refer to \code{\link[biomod2]{BIOMOD_Modeling}} to get more information about this metrics.
  If you give a \code{Fixed.thresh}, no optimisation will be done. Only the score for this threshold will be returned.

}
\value{
  A 1 row x 4 column \code{matrix} :
  \itemize{
    \item{\code{best.iter}:}{ the best score obtained for chosen statistic}
    \item{\code{cutoff}:}{ the associated cut-off used for transform fitted vector into binary}
    \item{\code{sensibility}:}{ the sensibility with this threshold}
    \item{\code{specificity}:}{ the specificity with this threshold}
  }
}


\author{ Damien Georges }

\seealso{ \code{\link[biomod2]{BIOMOD_Modeling}}, \code{\link[biomod2]{getStatOptimValue}}, \code{\link[biomod2]{calculate.stat}} }

\examples{
  a <- sample(c(0,1),100, replace=TRUE)
  
  # random drawing
  b <- runif(100,min=0,max=1000)
  Find.Optim.Stat(Stat='TSS',
                  Fit=b,
                  Obs=a,
                  Precision = 50,
                  Fixed.thresh = NULL)
  
  # biased drawing
  BiasedDrawing <- function(x, m1=300, sd1=200, m2=700, sd2=200){
    return(ifelse(x<0.5, rnorm(1,m1,sd1), rnorm(1,m2,sd2)))
  }
  
  c <- sapply(a,BiasedDrawing)
  
  Find.Optim.Stat(Stat='TSS',
                  Fit=c,
                  Obs=a,
                  Precision = 50,
                  Fixed.thresh = NULL)  
  
}

\keyword{ models, options, evaluate, evaluation }
