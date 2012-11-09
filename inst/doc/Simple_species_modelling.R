### R code from vignette source 'Simple_species_modelling.Rnw'
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
### code chunk number 2: loading_data
###################################################
# load the library
library(biomod2)

# load our species raster
# we consider only the presences of Myocastor coypus species
myResp.ras <- raster( system.file(
                        "external/species/Myocastor_coypus.img", 
                        package="biomod2") )

# extract the presences data

# the name
myRespName <- 'Myocastor'

# the XY coordinates of the presence
myRespXY <- xyFromCell(object=myResp.ras, 
                       cell=which(myResp.ras[]>0))

# and the presence data 
myResp <- extract(x=myResp.ras, y=myRespXY)

# load the environmental raster layers (could be .img, ArcGIS rasters or any supported format by the raster package)

# Environmental variables extracted from Worldclim (bio_3, bio_4, 
# bio_7, bio_11 & bio_12)
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
### code chunk number 3: formating_data
###################################################
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName,
                                     PA.nb.rep = 2,
                                     PA.nb.absences = 200,
                                     PA.strategy = 'random')


###################################################
### code chunk number 4: print_formating_data
###################################################
myBiomodData


###################################################
### code chunk number 5: plot_formating_data
###################################################
plot(myBiomodData)


###################################################
### code chunk number 6: modeling_options
###################################################
# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()


###################################################
### code chunk number 7: modeling
###################################################
# 3. Computing the models 

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


###################################################
### code chunk number 8: modeling_summary
###################################################
myBiomodModelOut 


###################################################
### code chunk number 9: modeling_model_evaluation
###################################################
# get all models evaluation                                     
myBiomodModelEval <- getModelsEvaluations(myBiomodModelOut)
                                     
# print the dimnames of this object
dimnames(myBiomodModelEval)
                                     
# let's print the TSS scores of Random Forest
myBiomodModelEval["TSS","Testing.data","RF",,]

# let's print the ROC scores of all selected models
myBiomodModelEval["ROC","Testing.data",,,]



###################################################
### code chunk number 10: modeling_variable_importance
###################################################
# print variable importances                                    
getModelsVarImport(myBiomodModelOut)


###################################################
### code chunk number 11: ensemble_modeling
###################################################
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


###################################################
### code chunk number 12: ensemble_modeling_outputs
###################################################
# print summary                     
myBiomodEM
                     
# get evaluation scores
getEMeval(myBiomodEM)


###################################################
### code chunk number 13: projection_curent
###################################################
# projection over the globe under current conditions
myBiomomodProj <- BIOMOD_Projection(
                         modeling.output = myBiomodModelOut,
                         new.env = myExpl,
                         proj.name = 'current',
                         selected.models = 'all',
                         binary.meth = 'ROC',
                         compress = 'xz',
                         clamping.mask = F)



###################################################
### code chunk number 14: projection_curent_plot
###################################################
# make some plots sub-selected by str.grep argument
plot(myBiomomodProj, str.grep = 'MARS')


###################################################
### code chunk number 15: projection_curent_getProj
###################################################
# if you want to make custom plots, you can also get the projected map
myCurrentProj <- getProjection(myBiomomodProj)
myCurrentProj


###################################################
### code chunk number 16: projection_future
###################################################
# load environmental variables for the future. 
myExpl2050 = stack( system.file( "external/climat/future/bio3.grd",
                                 package="biomod2"),
                    system.file( "external/climat/future/bio4.grd",
                                 package="biomod2"),
                    system.file( "external/climat/future/bio7.grd",
                                 package="biomod2"),
                    system.file( "external/climat/future/bio11.grd",
                                 package="biomod2"),
                    system.file( "external/climat/future/bio12.grd",
                                 package="biomod2"))

myBiomomodProj2050 <- BIOMOD_Projection(
                              modeling.output = myBiomodModelOut,
                              new.env = stack(myExpl2050),
                              proj.name = 't2050',
                              selected.models = 'all',
                              binary.meth = 'ROC',
                              compress = 'xz',
                              clamping.mask = T)
                              



###################################################
### code chunk number 17: projection_current_plot
###################################################
# make some plots, sub-selected by str.grep argument
plot(myBiomomodProj2050, str.grep = 'MARS')


###################################################
### code chunk number 18: EnsembleForecasting_future
###################################################
myBiomodEF <- BIOMOD_EnsembleForecasting( 
                      projection.output = myBiomomodProj2050,
                      EM.output = myBiomodEM )


###################################################
### code chunk number 19: EnsembleForecasting_loading_res
###################################################
load("Myocastor/proj_t2050/Myocastor_PA1_AllRun_EM.TSS")
Myocastor_PA1_AllRun_EM.TSS


###################################################
### code chunk number 20: EnsembleForecasting_plotting_res
###################################################
plot(Myocastor_PA1_AllRun_EM.TSS)


