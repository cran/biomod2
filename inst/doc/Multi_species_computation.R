### R code from vignette source 'Multi_species_computation.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: options
###################################################
options(prompt = " ", continue = "  ", width = 60, digits=4)
.CurFileName <- "biomod2_getting_started"
# .PrefixName <- strsplit(.CurFileName, "\\.")[[1]][1]
.PrefixName <- .CurFileName
.RversionName <- R.version.string
.PkgVersion <- packageDescription("biomod2")$Version


###################################################
### code chunk number 2: LoadSp_1
###################################################
# 1. loading species occurrences data
library(biomod2)

mySpeciesOcc <- read.csv( system.file( 
                          "external/species/species_occ.csv", 
                          package="biomod2"))
                            
head(mySpeciesOcc)


###################################################
### code chunk number 3: LoadEnv_1
###################################################
# 2. loading environmental data

# Environmental variables extracted from Worldclim (bio_3, bio_4, 
# bio_7, bio_11 & bio_12)
require(raster)
myExpl = stack( system.file( "external/climat/current/bio3.grd", 
                             package="biomod2"),
                system.file( "external/climat/current/bio4.grd", 
                             package="biomod2"), 
                system.file( "external/climat/current/bio7.grd", 
                             package="biomod2"),  
                system.file( "external/climat/current/bio11.grd", 
                             package="biomod2"), 
                system.file( "external/climat/current/bio12.grd", 
                             package="biomod2"))


###################################################
### code chunk number 4: Loop_1
###################################################

# define the species of interest
sp.names <- c("MelesMeles", "MyocastorCoypus")

# loop on species == applying the same functions to each species
for(sp.n in sp.names){
  
  cat('\n',sp.n,'modeling...')  
  ### definition of data 
  ## i.e keep only the column of our species
  myResp <- as.numeric(mySpeciesOcc[,sp.n])
  # get NAs id
  na.id <- which(is.na(myResp))
  # remove NAs to force the pseudo-absence extraction from background data 
  myResp <- myResp[-na.id]  ## presence-only data 
  
  myRespCoord = mySpeciesOcc[-na.id,c('x','y')]  ## coordinates of the presence-only data 
  
  myRespName = sp.n
  
  ### Initialisation
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       expl.var = myExpl,
                                       resp.xy = myRespCoord,
                                       resp.name = myRespName,
                                       PA.nb.rep = 2,
                                       PA.nb.absences = 10*sum(myResp==1, 
                                                               na.rm=TRUE),
                                       PA.strategy = 'random')
    
  
  ### Options definition
  myBiomodOption <- BIOMOD_ModelingOptions()
  
  ### Modelling 
  myBiomodModelOut <- BIOMOD_Modeling( 
                             myBiomodData, 
                             models = c('SRE','CTA','RF','MARS','FDA'), 
                             models.options = myBiomodOption, 
                             NbRunEval=1, 
                             DataSplit=80, 
                             Yweights=NULL, 
                             VarImport=3, 
                             models.eval.meth = c('TSS','ROC'),
                             SaveObj = TRUE,
                             rescal.all.models = TRUE)
  
  ### Building ensemble-models
  myBiomodEM <- BIOMOD_EnsembleModeling( 
                       modeling.output = myBiomodModelOut,
                       chosen.models = 'all',
                       eval.metric = c('TSS'),
                       eval.metric.quality.threshold = c(0.85),
                       prob.mean = T,
                       prob.cv = T,
                       prob.ci = T,
                       prob.ci.alpha = 0.05,
                       prob.median = T,
                       committee.averaging = T,
                       prob.mean.weight = T,
                       prob.mean.weight.decay = 'proportional' )
  
  ### Do projections on current variable
  myBiomomodProj <- BIOMOD_Projection(
                           modeling.output = myBiomodModelOut,
                           new.env = myExpl,
                           proj.name = 'current',
                           selected.models = 'all',
                           binary.meth= 'ROC',
                           compress = 'xz',
                           clamping.mask = F)
  
  ### Do ensemble-models projections on current variable
  myBiomodEF <- BIOMOD_EnsembleForecasting( 
                        projection.output = myBiomomodProj,
                        EM.output = myBiomodEM,
                        binary.meth = 'TSS',
                        total.consensus = TRUE)
}



###################################################
### code chunk number 5: alpha1
###################################################
# load the first species binary maps which will define the mask 
alphaMap <- get(load(paste(sp.names[1],"/proj_current/",
                           sp.names[1],"_TotalConsensus.bin.TSS",
                           sep="")))[[1]]

# free up the workspace to avoid memory saturation.  
rm(list=paste(sp.names[1],"_TotalConsensus.bin.TSS", sep=""))

# # add all other species map
for(sp.n in sp.names[-1]){
  # add layer
  alphaMap <- alphaMap + get(load(paste(sp.n,"/proj_current/",
                                        sp.n,"_TotalConsensus.bin.TSS",
                                        sep="")))[[1]]
  # free space
  rm(list=paste(sp.n,"_TotalConsensus.bin.TSS", sep=""))
}

# summary of created raster
alphaMap



###################################################
### code chunk number 6: alpha_2
###################################################
plot(alphaMap, main = expression( paste(alpha, "-diversity based on",
                                " TotalConsensus.bin.TSS outputs")))



###################################################
### code chunk number 7: Lapply_1
###################################################
MyBiomodSF <- function(sp.n){
  
  cat('\n',sp.n,'modeling...')  
  ### definition of data for this run
  ## i.e keep only the column of our species
  myResp <- as.numeric(mySpeciesOcc[,sp.n])
  # get NAs id
  na.id <- which(is.na(myResp))
  # remove NAs to enforce PA sampling to be done on explanatory rasters
  myResp <- myResp[-na.id]
  
  myRespCoord = mySpeciesOcc[-na.id,c('x','y')]
  
  myRespName = sp.n
  
  ### Initialisation
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       expl.var = myExpl,
                                       resp.xy = myRespCoord,
                                       resp.name = myRespName,
                                       PA.nb.rep = 2,
                                       PA.nb.absences = 10*sum(myResp==1,
                                                               na.rm=TRUE),
                                       PA.strategy = 'random')
    
  
  ### Options definition
  myBiomodOption <- BIOMOD_ModelingOptions()
  
  ### Modelling 
  myBiomodModelOut <- BIOMOD_Modeling( 
                             myBiomodData, 
                             models = c('SRE','CTA','RF','MARS','FDA'), 
                             models.options = myBiomodOption, 
                             NbRunEval=1, 
                             DataSplit=80, 
                             Yweights=NULL, 
                             VarImport=3, 
                             models.eval.meth = c('TSS','ROC'),
                             SaveObj = TRUE,
                             rescal.all.models = TRUE)
  
  ### Building ensemble-models
  myBiomodEM <- BIOMOD_EnsembleModeling( 
                       modeling.output = myBiomodModelOut,
                       chosen.models = 'all',
                       eval.metric = c('TSS'),
                       eval.metric.quality.threshold = c(0.85),
                       prob.mean = T,
                       prob.cv = T,
                       prob.ci = T,
                       prob.ci.alpha = 0.05,
                       prob.median = T,
                       committee.averaging = T,
                       prob.mean.weight = T,
                       prob.mean.weight.decay = 'proportional' )
  
  ### Do projections on current varaiable
  myBiomomodProj <- BIOMOD_Projection(
                           modeling.output = myBiomodModelOut,
                           new.env = myExpl,
                           proj.name = 'current',
                           selected.models = 'all',
                           binary.meth= 'ROC',
                           compress = 'xz',
                           clamping.mask = F)
  
  ### Do ensemble-models projections on current varaiable
  myBiomodEF <- BIOMOD_EnsembleForecasting( 
                        projection.output = myBiomomodProj,
                        EM.output = myBiomodEM,
                        binary.meth = 'TSS',
                        total.consensus = TRUE)
                        
}


###################################################
### code chunk number 8: Lapply_2
###################################################
myLapply_SFModelsOut <- lapply( sp.names, MyBiomodSF)


###################################################
### code chunk number 9: SnowFold_1 (eval = FALSE)
###################################################
## install.packages('snowfall', dependencies=TRUE)


###################################################
### code chunk number 10: SnowFold_2
###################################################
library(snowfall)


###################################################
### code chunk number 11: SnowFold_3 (eval = FALSE)
###################################################
## 
## ## Init snowfall
## library(snowfall)
## sfInit(parallel=TRUE, cpus=2 )  ## we select 2 CPUs. If you have 8 CPUs, put 8. 
## 
## ## Export packages
## sfLibrary('biomod2', character.only=TRUE)
## 
## ## Export variables
## sfExport('mySpeciesOcc')
## sfExport('myExpl')
## sfExport('sp.names')
## 
## # you may also use sfExportAll() to export all your workspace variables
## 
## ## Do the run
## mySFModelsOut <- sfLapply( sp.names, MyBiomodSF)
## 
## ## stop snowfall
## sfStop( nostop=FALSE )
## 


