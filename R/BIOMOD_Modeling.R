####################################################################################################
# BIOMOD_Modeling
# Damien.G
# feb 2012
####################################################################################################

# AIM :
#   Compute user's models selected and evaluate those.

# INPUT :
#   data <- a BIOMOD.formated.data object returned by BIOMOD_FormatingData() 
#   models <- vector of models names desired 
#   models.options <- a BIOMOD.models.options object returned by BIOMOD_ModelingOptions() 
#   NbRunEval <- Nb of Evaluation run 
#   DataSplit <- % of data used for models calibrations stuff
#   Yweights <- response points weights 
#   VarImport <- Nb of permutation done for variable Importances evaluation 
#   models.eval.meth <- vector of names of Models evaluation metrix 
#   SavePredictions <- keep or not array of predictions on hard disk (NOTE: Always TRUE for a posteriori EF)
#   KeepPredIndependent <- Not used
#   DoEnsembleForcasting <- Do or not EnsembleModeling.
#   SaveObj <- keep all results on hard disk or not (NOTE: strongly recommanded)

# OUTPUT : 
#   a BIOMOD.models.out (summary of all that were done) that will be given to others BIOMOD functions


# NOTE : 


####################################################################################################

'BIOMOD_Modeling' <- function( data, 
                               models = c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT'), 
                               models.options = NULL, 
                               NbRunEval=1, 
                               DataSplit=100, 
                               Yweights=NULL,
                               Prevalence=NULL,
                               VarImport=0, 
                               models.eval.meth = c('KAPPA','TSS','ROC'), 
#                                SavePredictions = TRUE, 
#                                KeepPredIndependent=FALSE, 
#                                DoEnsembleForcasting = TRUE,
                               SaveObj = TRUE,
                               rescal.all.models = TRUE,
                               do.full.models = TRUE){
  SavePredictions = TRUE
  # 0. loading required libraries =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  .Models.dependencies(silent=TRUE, models.options=models.options )

  # 1. args checking =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  args <- .Models.check.args(data, models, models.options, NbRunEval, DataSplit,
                             Yweights, VarImport, models.eval.meth, Prevalence)
  # updating Models arguments
  models <- args$models
  models.options <- args$models.options
  NbRunEval <- args$NbRunEval
  DataSplit <- args$DataSplit
  Yweights <- args$Yweights
  VarImport <- args$VarImport
  models.eval.meth <- args$models.eval.meth
  Prevalence <- args$Prevalence

  rm(args)
  models.out <- new('BIOMOD.models.out',
                    sp.name = data@sp.name,
                    expl.var.names = colnames(data@data.env.var),
                    has.evaluation.data = data@has.data.eval,
                    rescal.all.models = rescal.all.models)

#   #To keep track of Yweights state at origin (user's input)
#     if(NbRepPA!=0 && is.null(Yweights)) Yweights <- matrix(NA, nc=Biomod.material$NbSpecies, nr=nrow(DataBIOMOD))
   
  
  # 2. creating simulation directories =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #  
  # create the directories in which various objects will be stored (models, predictions and 
  # projection). The projection directories are created in the Projection() function.
  .Models.prepare.workdir(data@sp.name, models, SavePredictions)
  
  
  # 3. Saving Data and Model.option objects -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  if(SaveObj){
    # save Input Data
    save(data, file = paste(models.out@sp.name,"/.BIOMOD_DATA/formated.input.data",sep=""), compress=ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz'))
    models.out@formated.input.data@inMemory <- FALSE
    models.out@formated.input.data@link <- paste(models.out@sp.name,"/.BIOMOD_DATA/formated.input.data",sep="")
    # save Model Options
    save(models.options, file = paste(models.out@sp.name,"/.BIOMOD_DATA/models.options",sep=""), compress=ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz'))
    models.out@models.options@inMemory <- FALSE
    models.out@models.options@link <- paste(models.out@sp.name,"/.BIOMOD_DATA/models.options",sep="")

  }
  

  # 3. rearanging data and determining calib and eval data -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

  mod.prep.dat <- .Models.prepare.data(data, NbRunEval, DataSplit, Yweights, Prevalence, do.full.models)
  rm(data)
  
  # keeping calibLines
#   if(NbRunEval > 0){
    calib.lines <- mod.prep.dat[[1]]$calibLines
    if(length(mod.prep.dat) > 1){
      for(pa in 2:length(mod.prep.dat)){
        calib.lines <- abind(calib.lines, mod.prep.dat[[pa]]$calibLines, along=3)
      }
    }
    # save calib.lines
    save(calib.lines, file = paste(models.out@sp.name,"/.BIOMOD_DATA/calib.lines",sep=""), compress=ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz'))
    models.out@calib.lines@inMemory <- FALSE
    models.out@calib.lines@link <- paste(models.out@sp.name,"/.BIOMOD_DATA/calib.lines",sep="")
#     models.out@calib.lines@val <- calib.lines
    rm(calib.lines)
#   } 

  # 4. Print modelling summary in console -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  .Models.print.modeling.summary(mod.prep.dat, models)
  
  # 5. Runing models -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  # loop on PA
  modeling.out <- lapply(mod.prep.dat,.Biomod.Models.loop,
                      Model = models,
                      Options = models.options,
                      VarImport = VarImport, 
                      mod.eval.method = models.eval.meth,
#                       do.EF = DoEnsembleForcasting,
                      SavePred = SavePredictions,
                      rescal.models = rescal.all.models
                      )
  # put outputs in good format and save those

  models.out@models.computed <- .transform.outputs(modeling.out, out='models.run')
  models.out@models.failed <- .transform.outputs(modeling.out, out='calib.failure')
  
  if(SaveObj){
    # save model evaluation
    models.evaluation <- .transform.outputs(modeling.out, out='evaluation')
    save(models.evaluation, file = paste(models.out@sp.name,"/.BIOMOD_DATA/models.evaluation",sep=""), compress=ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz'))
    models.out@models.evaluation@inMemory <- TRUE
    models.out@models.evaluation@link <- paste(models.out@sp.name,"/.BIOMOD_DATA/models.evaluation",sep="")
    models.out@models.evaluation@val <- models.evaluation
    rm(models.evaluation)
    
    # save model variables importances
    if(VarImport > 0 ){
      variables.importances <- .transform.outputs(modeling.out, out='var.import')
      ## trick to put appropriate dimnames
      vi.dim.names <- dimnames(variables.importances)
      vi.dim.names[[1]] <- models.out@expl.var.names
      dimnames(variables.importances) <- vi.dim.names
      rm('vi.dim.names')
    
      save(variables.importances, file = paste(models.out@sp.name,"/.BIOMOD_DATA/variables.importances",sep=""), compress=ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz'))
      models.out@variables.importances@inMemory <- TRUE
      models.out@variables.importances@link <- paste(models.out@sp.name,"/.BIOMOD_DATA/variables.importances",sep="")
      models.out@variables.importances@val <- variables.importances
      rm(variables.importances)
    }

    # save model predictions
    models.prediction <- .transform.outputs(modeling.out, out='prediction')
    save(models.prediction, file = paste(models.out@sp.name,"/.BIOMOD_DATA/models.prediction",sep=""), compress=ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz'))
    models.out@models.prediction@inMemory <- FALSE
    models.out@models.prediction@link <- paste(models.out@sp.name,"/.BIOMOD_DATA/models.prediction",sep="")
#     models.out@models.prediction@val <- .transform.outputs(modeling.out, out='prediction')
    rm(models.prediction)
    
    # save evaluation model predictions
    models.prediction.eval <- .transform.outputs(modeling.out, out='prediction.eval')
    save(models.prediction.eval, file = paste(models.out@sp.name,"/.BIOMOD_DATA/models.prediction.eval",sep=""), compress=ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz'))
    models.out@models.prediction.eval@inMemory <- FALSE
    models.out@models.prediction.eval@link <- paste(models.out@sp.name,"/.BIOMOD_DATA/models.prediction.eval",sep="")
#     models.out@models.prediction@val <- .transform.outputs(modeling.out, out='prediction')
    rm(models.prediction.eval)

  }
  
  # removing MAXENT tmp dir
  if('MAXENT' %in% models){
    .Delete.Maxent.WorkDir(species.name=models.out@sp.name)
  }
  
  rm(modeling.out)
  
  # create models id
  eval(parse(text=paste(models.out@sp.name,'.models.out <- models.out', sep='')))
  eval(parse(text=paste("save(",models.out@sp.name,".models.out,file= '",
                         models.out@sp.name,"/",models.out@sp.name,".models.out')"
                        ,sep="")))

  cat("\n-=-=-=- Done -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")
  return(models.out)
}

# -=-=-=- Several hidden functions -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #


.Models.dependencies <- function(silent=TRUE, models.options = NULL){
  # Loading all required libraries
  cat('\n\nLoading required library...')
  require(nnet, quietly=silent)
  require(rpart, quietly=silent)
  require(MASS, quietly=silent)
  require(gbm, quietly=silent)
  require(mda, quietly=silent)
  require(randomForest, quietly=silent)
  
  if(!is.null(models.options)){
    if(grepl('mgcv', models.options@GAM$algo)){
      if("package:gam" %in% search() ) detach(package:gam)
      require(mgcv, quietly=silent)
    } else{
      if("package:mgcv" %in% search() ) detach(package:mgcv)
      require(gam, quietly=silent)
    }
  } else {
    if('mgcv' %in% rownames(installed.packages())){
      if("package:gam" %in% search() ) detach(package:gam)
      require(mgcv, quietly=silent)
    } else{
      if("package:mgcv" %in% search() ) detach(package:mgcv)
      require(gam, quietly=silent)
    }    
  }

  require(abind, quietly=silent)
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.Models.check.args <- function(data, models, models.options, NbRunEval, DataSplit,
                               Yweights, VarImport, models.eval.meth, Prevalence){
  cat('\n\nChecking Models arguments...\n')
  # data checking
  if( !(class( data ) %in% c("BIOMOD.formated.data",
                             "BIOMOD.formated.data.PA",
                             "BIOMOD.formated.data.indep",
                             "BIOMOD.formated.data.PA.indep") )){
    stop( "data argument mut be a 'BIOMOD.formated.data' (obtain by runing Initial.State function) ")
  }
  
  # models checking
  if( !is.character( models ) ){
    stop("models argument must be a 'character' vector")
  }
  
  models <- unique(models)
  
  if(sum(models %in% c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT')) != length(models)){
    stop(paste(models[which( (models %in% c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT'))
                             == FALSE) ]," is not a availabe model !",sep=""))
  }
  
  # models.options checking ( peut etre permetre l'utilisation de liste de params )
  if( !is.null(models.options) && class(models.options) != "BIOMOD.Model.Options" ){
    stop("models.options argument must be a 'BIOMOD.Model.Options.object' (obtain by runing ... ) ")
  }
  
  if( is.null(models.options)){
    warning("Models will run with 'defaults' parameters", immediate.=T)
    # create a default models.options object
    models.options = new("BIOMOD.Model.Options")
  }
  
  # MAXENT specific checking
  if("MAXENT" %in% models){
    if(!file.exists(file.path(models.options@MAXENT$path_to_maxent.jar ,"maxent.jar")) ){
      models = models[-which(models=='MAXENT')]
      warning("The maxent.jar file is missing. You need to download this file (http://www.cs.princeton.edu/~schapire/maxent) and put the maxent.jar file in your working directory -> MAXENT was switched off")
    }
    if(!.check.java.installed()){
      models = models[-which(models=='MAXENT')]
    } else{
      if(nrow(data@coord)==1){ # no coordinates
        warning("You must give XY coordinates if you want to run MAXENT -> MAXENT was switched off")
        models = models[-which(models=='MAXENT')]
      }
    } 
  }
  
  # NbRunEval checking (permetre un nb different par espece?)
  if( !is.numeric(NbRunEval) || NbRunEval <= 0 ){
    stop("NbRunEval argument mus be a non null positive 'numeric'")
  }
  
  # DataSplit checking
  if( !is.numeric(DataSplit) || DataSplit < 0 || DataSplit > 100 ){
    stop("DataSplit argument must be a 0-100 'numeric'")
  }
  
  if(DataSplit < 50){
    warning("You choose to allocate more data to evaluation than to calibration of your model 
            (DataSplit<50)\nMake sure you really wanted to do that. \n", immediate.=T)
  }
  
#   # EM weight checking
#   if(!is.null(EM.weight))
#     if(!any(EM.weight==c("Roc","TSS","Kappa")))
#       stop("The 'EM.weight' parameter must be one of the following: NULL, 'Roc', 'TSS' or 'Kappa'.\n")
#       
  # Check that the weight matrix was entered correctly  
  if(!is.null(Yweights)){
     if(!is.numeric(Yweights))
        stop("Yweights must be a numeric vector")
     if(length(Yweights) != length(data@data.species)) 
       stop("The number of 'Weight' does not match with the input calibration data. 
            Simulation cannot proceed.")
  }
  
  # Defining evaluation runs.
  if(NbRunEval <= 0){
      DataSplit <- 100
      if(!(class(data) %in% c("BIOMOD.formated.data.indep", "BIOMOD.formated.data.PA.indep") )){
        warning("The models will be evaluated on the calibration data only (NbRunEval=0 and no
                independent data) \n\t it could lead to over-optimistic predictive performances.\n",
                immediate.=T)
      }
  }
  if(DataSplit==100) NbRunEval <- 0
  
  # Switch SRE and MARS off if one of the variables is a non-numeric.
  if(sum(models %in% c('MARS', 'SRE')) > 0 ){
    if(sum(apply(data@data.env.var,2,is.factor)) > 0){
      modId <- which((models %in% c('MARS', 'SRE') ) == TRUE)
      warning(paste("Because some environmental variables are factors, ",models[modId],
                    " have been switched off! \n", sep=""))
      models <- models[-modId] 
    }
  }
  
  
  # Models evaluation method checking
  models.eval.meth <- unique(models.eval.meth)
  
  if(sum(models.eval.meth %in% c('FAR','SR','HSS','ORSS','TSS','KAPPA','ACCURACY','BIAS',
                              'POD','PODFD','CSI','ETS','HK','ROC')) != length(models.eval.meth)){
    stop(paste(models.eval.meth[which( (models.eval.meth %in% c('FAR','SR','HSS','ORSS','TSS',
                                                                'KAPPA','ACCURACY','BIAS', 'POD',
                                                                'PODFD','CSI', 'ETS','HK','ROC')) 
                                       == FALSE) ]," is not a availabe models evaluation metric !",sep=""))
  }
  
  # Prevalence checking
  if(!is.null(Prevalence)){
    if(!is.numeric(Prevalence) | Prevalence>=1 | Prevalence <=0){
      stop("Prevalence must be a 0-1 numeric")
    }    
  }

#   cat('\nChecking done!\n')
  return(list(models = models,
              models.options = models.options,
              NbRunEval = NbRunEval,
              DataSplit = DataSplit,
              Yweights = Yweights,
              VarImport = VarImport,
              models.eval.meth = models.eval.meth,
              Prevalence = Prevalence))
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.Models.prepare.workdir <- function(sp.name, model, SavePredictions=FALSE){
  cat("\nCreating suitable Workdir...\n")
  dir.create(paste(getwd(),"/",sp.name,sep=""), showWarnings=FALSE)
  dir.create(paste(getwd(),"/",sp.name,"/.BIOMOD_DATA",sep=""), showWarnings=FALSE)
  dir.create(paste(getwd(),"/",sp.name, "/models", sep=""), showWarnings=FALSE)
#   if(SavePredictions){
#     dir.create(paste(getwd(),"/",sp.name, "/pred", sep=""), showWarnings=FALSE)
#   }
  if(sum(model %in% c('MARS', 'FDA', 'ANN')) > 0 ){
    dir.create(paste(getwd(),"/",sp.name, "/models/rescaling_models", sep=""), showWarnings=FALSE)
  }
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
.SampleMat <- function(data, dataSplit, nbRun = 1){
  # return a matrix with nbRun columns of boolean (T: calib, F= eval)
  # data is a 0,1 vector
  pres <- which(data==1)
  abs <- (1:length(data))[-pres]
  
  nbPresEval <- round(length(pres) * dataSplit/100)
  nbAbsEval <- round(length(abs) * dataSplit/100)
  
  mat.out <- matrix(FALSE, nrow=length(data), ncol=nbRun)
  colnames(mat.out) <- paste('_RUN',1:nbRun, sep='')
  
  for (i in 1:ncol(mat.out)){
    mat.out[sample(pres, nbPresEval),i] <- TRUE
    mat.out[sample(abs, nbAbsEval),i] <- TRUE
  }
  
  return(mat.out)
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.Models.print.modeling.summary <- function( mod.prep.dat, models){
  cat("\n\n-=-=-=- ", unlist(strsplit(mod.prep.dat[[1]]$name,'_'))[1], " Modeling Summary -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n")

  cat(ncol(mod.prep.dat[[1]]$dataBM)-1, " environmental variables (", colnames(mod.prep.dat[[1]]$dataBM)[-1], ")\n")

  cat("Number of evaluation repetitions :" , ncol(mod.prep.dat[[1]]$calibLines), "\n")

  cat("Models selected :", models, "\n")

  cat("Total number of model runs :",ncol(mod.prep.dat[[1]]$calibLines) * length(models) * length(mod.prep.dat),"\n")
  
  cat("-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")
}
  
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  
.Models.check.EF.args <- function(models, models.eval.meth, models.options){
  # the models selected for doing EF
  if(models.options@EF$models.selected == 'all'){
    EF.algo <- models[which(models != 'SRE')] # remove SRE technique if it was selected (SRE cannot be used for ensemble forecast)
  } else {
    EF.algo <- models[models %in% models.options@EF$models.selected]
  }
  if(length(EF.algo)==0) stop('No models available selected for Ensemble forcastiong stuff')
  # the weight methods
  if(models.options@EF$weight.method == 'all'){
    EF.weight <- models.eval.meth
  } else {
    EF.weight <- models.eval.meth[models.eval.meth %in% models.options@EF$weight.method]
    cat('\n***', EF.weight)
  }
  if(length(EF.weight)==0) stop('No weighting method available selected for Ensemble forcastiong stuff')
  return(list(EF.algo = EF.algo,
            EF.weight = EF.weight))
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
.check.java.installed <- function(){
  
  if(.Platform$OS.type == "unix"){
    java.test <- try( expr = eval(system("command -v java", intern=TRUE )) , silent=TRUE)
  } else if(.Platform$OS.type == "windows"){
    java.test <- try( expr = eval(system( "java", intern=TRUE )), silent=TRUE )
  } else java.test <- ""
  
  if(!is.null(attr(java.test,"class"))){
    cat("\n! java software seems not be corectly installed\n  > MAXENT modelling was switched off!")
    return(FALSE)
  } else{ return(TRUE) }
  
}

