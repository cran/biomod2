\name{BIOMOD_EnsembleForecasting}
\alias{BIOMOD_EnsembleForecasting}

\title{ Ensemble projections of species over space and time }
\description{
This function use projections of \sQuote{individual models} and ensemble models from code{\link[biomod2]{BIOMOD_EnsembleModeling}} to build an ensemble of species' projections over space and time.
}
\usage{
BIOMOD_EnsembleForecasting( projection.output,
                            EM.output,
                            total.consensus = FALSE,
                            binary.meth = NULL,
                            filtered.meth = NULL,
                            ...)
}


\arguments{
  \item{projection.output}{a \code{"\link[=BIOMOD.projection.out-class]{BIOMOD.projection.out}"} returned by \code{\link[biomod2]{BIOMOD_Projection}} }
  \item{EM.output}{ a \code{"\link[=BIOMOD.EnsembleModeling.out-class]{BIOMOD.EnsembleModeling.out}"} returned by \code{\link[biomod2]{BIOMOD_EnsembleModeling}}}
  \item{total.consensus}{Deprecated, construct ensemble-models setting \code{em.by} argument to \code{'all'} to get total consensus-models}
  \item{binary.meth}{vector specifying the names of evaluation metrics and associated thresholds to transform the probabilities of presence into presence and absence (binary transformation).  }
  \item{filtered.meth}{vector specifying the names of evaluation metrics and associated thresholds to transform into 0 the probabilities of presence lower than the thresholds.}
  \item{\ldots}{ further arguments (see details)}
}

\details{
This function requires to have successfully run \pkg{biomod2} modeling, ensemble-modeling and projection steps. Ensemble projections will be created in respect to \code{projection.output} projections, which are combined following \code{EM.output} ensemble-modeling rules.

The \sQuote{total.consensus} projection is basically the mean of all projections (for having only one output).  
} % end of detail 

\value{
Nothing returned but specific \sQuote{projection files} are saved on the hard drive projection folder. This files are either an \code{array} or a \code{RasterStack} depending the original projections data type.
Load these created files to plot and work with them.

}% end of value


\author{ 
Wilfried Thuiller, Damien Georges, Robin Engler
}

\seealso{
\code{\link[biomod2]{BIOMOD_EnsembleModeling}}, \code{\link[biomod2]{BIOMOD_Projection}}
}

\examples{
# 0. Load data & Selecting Data

# species occurrences
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
                                    package="biomod2"), row.names = 1)
head(DataSpecies)

# the name of studied species
myRespName <- 'GuloGulo'

# the presence/absences data for our species 
myResp <- as.numeric(DataSpecies[,myRespName])

# the XY coordinates of species data
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]


# Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
myExpl = stack( system.file( "external/bioclim/current/bio3.grd", 
                             package="biomod2"),
                system.file( "external/bioclim/current/bio4.grd", 
                             package="biomod2"), 
                system.file( "external/bioclim/current/bio7.grd", 
                             package="biomod2"),  
                system.file( "external/bioclim/current/bio11.grd", 
                             package="biomod2"), 
                system.file( "external/bioclim/current/bio12.grd", 
                             package="biomod2"))

# 1. Formatting Data
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)
                                                                     
# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

# 3. Running the models 
myBiomodModelOut <- BIOMOD_Modeling( myBiomodData, 
                                       models = c('SRE','RF'), 
                                       models.options = myBiomodOption, 
                                       NbRunEval=2, 
                                       DataSplit=80, 
                                       Yweights=NULL, 
                                       VarImport=0, 
                                       models.eval.meth = c('TSS'),
                                       SaveObj = TRUE,
                                       rescal.all.models = FALSE,
                                       do.full.models = FALSE)
                                       
# 4. Creating the ensemble models 
myBiomodEM <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut,
                                       chosen.models = grep('_RF', getModelsBuiltModels(myBiomodModelOut), value=TRUE),
                                       em.by = 'algo',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.7),
                                       prob.mean = FALSE,
                                       prob.cv = FALSE,
                                       prob.ci = FALSE,
                                       prob.ci.alpha = 0.05,
                                       prob.median = FALSE,
                                       committee.averaging = FALSE,
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional' )
                                       
# 5. Individual models projections on current environmental conditions
myBiomodProjection <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                        new.env = myExpl,
                                        proj.name = 'current',
                                        selected.models = grep('_RF', getModelsBuiltModels(myBiomodModelOut), value=TRUE),
                                        compress = 'NULL',
                                        build.clamping.mask = FALSE)
                                          

# 4. Creating the ensemble projections
BIOMOD_EnsembleForecasting( projection.output = myBiomodProjection,
                            EM.output = myBiomodEM)
}

\keyword{ models }
