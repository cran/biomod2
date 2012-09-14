\name{BIOMOD_Projection}
\alias{BIOMOD_Projection}

\title{ Project the calibrated models within \pkg{biomod2} into new space or time}
\description{
  For all the models currently implemented, \pkg{biomod2} is able to project potential distributions of species in other areas, other resolutions or other time scales.
}
\usage{
BIOMOD_Projection(modeling.output,
                  new.env,
                  proj.name,
                  xy.new.env = NULL,
                  selected.models = 'all',
                  binary.meth = NULL,
                  filtered.meth = NULL,
                  compress = 'xz',
                  clamping.mask = TRUE,
                  ...)
}

\arguments{

\item{modeling.output}{  \code{"\link[=BIOMOD.models.out-class]{BIOMOD.models.out}"} object produced by a \code{\link[biomod2]{BIOMOD_Modeling}} run}

\item{new.env}{ A set of explanatory variables onto which models will be projected . It could be a \code{data.frame}, a \code{matrix}, or a \code{rasterStack} object. Make sure the column names (\code{data.frame} or \code{matrix}) or layerNames (\code{rasterStack}) perfectly match with the names of variables used to build the models in the previous steps.}

\item{proj.name}{a character defining the projection name (a new folder will be created with this name)}

\item{xy.new.env}{optional coordinates of new.env data. Ignored if new.env is a \code{rasterStack}}

\item{selected.models}{'all' when all models have to be used to render projections or a subset vector of modeling.output models computed (accessing with the slot @models.computed of your \code{"\link[=BIOMOD.models.out-class]{BIOMOD.models.out}"} object)}

\item{binary.meth}{a vector of a subset of models evaluation method computed before (see \code{\link[biomod2]{BIOMOD_Modeling}}). If NULL then no binary transformation computed, else the given binary techniques will be used to transform the projection into 0/1 data. }

\item{filtered.meth}{a vector of a subset of models evaluation method computed before (see \code{\link[biomod2]{BIOMOD_Modeling}}). if NULL then no filtering transformation computed, else the given binary techniques will be used to transform the projection by settting to 0 the probability values below a specific threshold.}

\item{compress}{compression format of objects stored on your hard drive. May be one of \sQuote{xz}, \sQuote{gzip} or \code{NULL}}

\item{clamping.mask}{ if TRUE, a clamping mask will be saved on hard drive different (see details)}

\item{...}{ Additional arguments }
}


\details{
Projections are done for all selected models, that means (by default) for all evaluation run, and pseudo absences selections if applicable. This projections may be used later to compute ensemble forecasting.

If \code{clamping.mask} is set to \code{TRUE} a file (same type than \code{new.env} arg) will be saved in your projection folder. This mask will identifies locations where predictions are uncertain because the values of the variables are outside the range used for calibrating the models. The \sQuote{clamping.mask} values correrspond to the number of variables that are out of their calibrating/training range. (see vignette for more details)

\ldots may be :
\itemize{
  \item{\code{do.stack}: }{logical. return and save projection in a lone \code{RasterStack} if it's possible.}
}

}

\value{
Returns the projections for all selected model (\code{"\link[=BIOMOD.projection.out-class]{BIOMOD.projection.out}"} object), and stored in the hard drive on the specific directory names by the name of the projection. The data is a 4-dimensions array (see ...) if new.env is a \code{matrix} or a \code{data.frame}. It is a rasterStack if new.env is a \code{rasterStack} and or several rasterLayers if the \code{rasterStack} is too large.

A new folder is also created on your hard drive. This folder contains the created projection object (basic one and binary and filtered ones if selected). The object are loaded with the \code{\link{load}} function. The loaded object can be then plotted and analyzed.

}

\author{
Wilfried Thuiller, Damien Georges
}

\seealso{
\code{\link{BIOMOD_Modeling}, \link{BIOMOD_FormatingData}},  \code{\link{BIOMOD_ModelingOptions}}
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
                                     PA.nb.rep = 1,
                                     PA.nb.absences = 200,
                                     PA.strategy = 'random')
                                                                     
# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

# 3. Doing Modelisation

myBiomomodModelOut <- BIOMOD_Modeling( myBiomodData, 
                                       models = c('CTA','RF'), 
                                       models.options = myBiomodOption, 
                                       NbRunEval=1, 
                                       DataSplit=80, 
                                       Yweights=NULL, 
                                       VarImport=3, 
                                       models.eval.meth = c('TSS','ROC'),
                                       SaveObj = TRUE )
                                       
                                       
# 4.1 Projection on current environemental conditions

myBiomomodProjection <- BIOMOD_Projection(modeling.output = myBiomomodModelOut,
                                          new.env = myExpl,
                                          proj.name = 'current',
                                          selected.models = 'all',
                                          binary.meth = 'ROC',
                                          filtered.meth = 'TSS',
                                          compress = 'xz',
                                          clamping.mask = TRUE)
                                          
# 4.2 Projection on future environemental conditions

myExpl2050 = raster::stack(system.file("external/climat/future/bio3.grd",package="biomod2"),
                           system.file("external/climat/future/bio4.grd",package="biomod2"),
                           system.file("external/climat/future/bio7.grd",package="biomod2"),
                           system.file("external/climat/future/bio11.grd",package="biomod2"),
                           system.file("external/climat/future/bio12.grd",package="biomod2"))

myBiomomodProjection2050 <- BIOMOD_Projection(modeling.output = myBiomomodModelOut,
                                              new.env = myExpl2050,
                                              proj.name = 't2050',
                                              selected.models = 'all',
                                              binary.meth = 'ROC',
                                              filtered.meth = 'TSS',
                                              compress = 'xz',
                                              clamping.mask = TRUE)

# print summary and plot projections
myBiomomodProjection2050
plot(myBiomomodProjection2050)
                                       
}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ tree }

