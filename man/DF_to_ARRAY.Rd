\name{DF_to_ARRAY}
\alias{DF_to_ARRAY}

\title{Convert a biomod2 data.frame into array }

\description{ DF_to_ARRAY is an internal \pkg{biomod2} function that can be useful to help users to transform a standard \pkg{biomod2} \code{data.frame} output into the standard \code{array} one.
}

\usage{
DF_to_ARRAY(df)
}



\arguments{
  \item{df}{  a standard \pkg{biomod2} \code{data.frame} output }
}


\details{
This function can be useful when you want to re-convert \pkg{biomod2} \code{data.frame} (e.g. projections) into their initial \code{array} format (see \code{\link[biomod2]{BIOMOD_Projection}} for further details) 
}


\value{
A standard \pkg{biomod2} \code{array} output.
}


\author{ Damien Georges }

\seealso{ \code{\link[biomod2]{BIOMOD_Projection}} }
\examples{
# 0. Load data & Selecting Data

# species occurrences
species_occ <- read.csv(system.file("external/species/species_occ.csv",package="biomod2"))

# we consider only presences of MyocastorCoypus species
myRespName <- 'MyocastorCoypus'
myRespCoord <- species_occ[,c('x','y')]
myResp <- as.numeric(species_occ[,myRespName])

# Environemental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
myExpl = raster::stack(system.file("external/climat/current/bio3.grd",package="biomod2"),
                       system.file("external/climat/current/bio4.grd",package="biomod2"),
                       system.file("external/climat/current/bio7.grd",package="biomod2"),
                       system.file("external/climat/current/bio11.grd",package="biomod2"),
                       system.file("external/climat/current/bio12.grd",package="biomod2"))

# Keep only points where we have info                       
myExpl <- extract(myExpl, myRespCoord)

# 1. Formating Data
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespCoord,
                                     resp.name = myRespName,
                                     PA.nb.rep = 1,
                                     PA.nb.absences = 100,
                                     PA.strategy = 'random')
                                                                     
# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

# 3. Doing Modelisation

myBiomomodModelOut <- BIOMOD_Modeling( myBiomodData, 
                                       models = c('SRE','RF'), 
                                       models.options = myBiomodOption, 
                                       NbRunEval=1, 
                                       DataSplit=80, 
                                       Yweights=NULL, 
                                       VarImport=3, 
                                       models.eval.meth = c('TSS'),
                                       SaveObj = TRUE )
                                       
                                       
# 4 Projection on current environemental conditions

myBiomomodProjection <- BIOMOD_Projection(modeling.output = myBiomomodModelOut,
                                          new.env = data.frame(myExpl),
                                          proj.name = 'current',
                                          selected.models = 'all')
                                          


# 5. Get projection under data.frame format
myProjDF <- getProjection(myBiomomodProjection, as.data.frame=TRUE)
class(myProjDF)
dim(myProjDF)
dimnames(myProjDF)

# 6. Transform data.frame into array
myProjArray <- DF_to_ARRAY(myProjDF)
class(myProjArray)
dim(myProjArray)
dimnames(myProjArray)

}

\keyword{ models, formula, options }
