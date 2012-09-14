# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# Definition of models build with biomod2 to make it easy to access
# plot, .predict...
# Damien G. - 20/11/2012
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# NOTE -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# This script will induce strong changes in biomod2. The version
# of biomod2 using objects defined there will be biomod2_2.x.x
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #


# Formal Class =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('biomod2_model',
         representation(model_name = 'character',
                        model_class = 'character',
                        model_options = 'list',
                        model = 'ANY',
                        scaling_model = 'ANY',
                        resp_name = 'character',
                        expl_var_names = 'character',
                        model_evaluation = 'matrix',
                        model_variables_importance = 'matrix'),
         prototype(model_name = 'mySpecies_DataSet_RunName_myModelClass',
                   model_class = 'myModelClass',
                   model_options = list(),
                   model = list(),
                   scaling_model = list(),
                   resp_name = 'mySpecies',
                   expl_var_names = 'myRespVar',
                   model_evaluation = matrix(),
                   model_variables_importance = matrix()),
         validity = function(object){
           
           # check that scaler is a glm if it is defined
           if(length(object@scaling_model)) 
             if(sum(! ( c("glm", "lm") %in% class(object@scaling_model) ) ) > 0) 
               return(FALSE)
           
           return(TRUE) 
           } )

setMethod('show', signature('biomod2_model'),
          function(object){
            .bmCat("'biomod2_model'")
            cat("\n\t model name :", object@model_name, fill=.Options$width)
            cat("\n\t model class :", object@model_class, fill=.Options$width)
            cat("\n\t This model", ifelse(length(object@scaling_model), "has", "does'nt have"),"its own scaler", fill=.Options$width)
            
            cat("\n")
            cat("\n\t response modelled :", object@resp_name, fill=.Options$width)
            cat("\n\t explanatory variables used:", object@expl_var_names, fill=.Options$width)
            
            cat("\n")
            cat("\n\t NOTE : ")
            cat("\n\t\t You can access 'formal' model with getFormalModel function")
            cat(ifelse(length(object@scaling_model), "\n\t\t You can access scaling model with getScalingModel function\n", "\n"))

            .bmCat()
          })

# Models getters -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
setGeneric( "getFormalModel", 
            def = function(object){
              standardGeneric( "getFormalModel" )
            } )

setMethod('getFormalModel', signature('biomod2_model'),
          function(object){
            return(object@model)
          })

setGeneric( "getScalingModel", 
            def = function(object){
              standardGeneric( "getScalingModel" )
            } )

setMethod('getScalingModel', signature('biomod2_model'),
          function(object){
            return(object@scaling_model)
          })

# ANN Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

setClass('ANN_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype(model_class = 'ANN'),
         validity = function(object){
           # check model class
           if(sum(! ( c("nnet.formula", "nnet") %in% class(object@model) ) ) > 0) return(FALSE)
           return(TRUE)
           })

setMethod('predict', signature(object = 'ANN_biomod2_model'),
          function(object, newdata, ...){
            
            args <- list(...)
            
            if( ! ("package:nnet" %in% search()) ){ require(nnet,quietly=TRUE) }
            
            if(inherits(newdata, 'Raster')){            
              return(.predict.ANN_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
               return(.predict.ANN_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }
            
          })

.predict.ANN_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000
  
  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  
  set.seed(555)
  proj <- predict(newdata, getFormalModel(object), type="raw")
  
  if(length(getScalingModel(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = getScalingModel(object), Prev = 0.5 , dat = proj)
  }
  
  if(on_0_1000) proj <- round(proj*1000)
  
  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    writeRaster(proj, filename=filename, overwrite=overwrite)
    proj <- raster(filename)
  }

  return(proj)
}

.predict.ANN_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  
  set.seed(555)
  proj <- as.numeric( predict(getFormalModel(object), newdata, type="raw") )
  
  if(length(getScalingModel(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = getScalingModel(object), Prev = 0.5 , dat = proj)
  }
  
  if(on_0_1000) proj <- round(proj*1000)
  
  return(proj)
}





# CTA Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('CTA_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype(model_class = 'CTA'),
         validity = function(object){
           # check model class
           if(sum(! ( c("rpart") %in% class(object@model) ) ) > 0) return(FALSE)
           return(TRUE)
         })

setMethod('predict', signature(object = 'CTA_biomod2_model'),
          function(object, newdata, ...){
            
            args <- list(...)
            
            if( ! ("package:rpart" %in% search()) ){ require(rpart,quietly=TRUE) }
            
            if(inherits(newdata, 'Raster')){            
              return(.predict.CTA_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.CTA_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }
            
          })

.predict.CTA_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000
  
  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  
  set.seed(123) # to be able to refind our trees MAY BE BAD
  proj <- predict(newdata, model=getFormalModel(object), type='prob', index=2)
  
  if(length(getScalingModel(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = getScalingModel(object), Prev = 0.5 , dat = proj)
  }
  
  if(on_0_1000) proj <- round(proj*1000)
  
  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    writeRaster(proj, filename=filename, overwrite=overwrite)
    proj <- raster(filename)
  }
  
  return(proj)
}

.predict.CTA_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  
  set.seed(123)
  proj <- as.numeric(predict(getFormalModel(object), as.data.frame(newdata),type="prob")[,2])
  
  if(length(getScalingModel(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = getScalingModel(object), Prev = 0.5 , dat = proj)
  }
  
  if(on_0_1000) proj <- round(proj*1000)
  
  return(proj)
}





# FDA Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('FDA_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype(model_class = 'FDA'),
         validity = function(object){
           # check model class
           if(sum(! ( c("fda") %in% class(object@model) ) ) > 0) return(FALSE)
           return(TRUE)
         })

setMethod('predict', signature(object = 'FDA_biomod2_model'),
          function(object, newdata, ...){
            
            args <- list(...)
            
            if( ! ("package:mda" %in% search()) ){ require(mda,quietly=TRUE) }
            
            if(inherits(newdata, 'Raster')){            
              return(.predict.FDA_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.FDA_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }
            
          })

.predict.FDA_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000
  
  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  
  proj <- predict(newdata, model=getFormalModel(object), type='posterior', index=2)
  
  if(length(getScalingModel(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = getScalingModel(object), Prev = 0.5 , dat = proj)
  }
  
  if(on_0_1000) proj <- round(proj*1000)
  
  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    writeRaster(proj, filename=filename, overwrite=overwrite)
    proj <- raster(filename)
  }
  
  return(proj)
}

.predict.FDA_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  
  proj <- as.numeric(predict(getFormalModel(object), as.data.frame(newdata),type = "posterior")[,2])

  if(length(getScalingModel(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = getScalingModel(object), Prev = 0.5 , dat = proj)
  }
  
  if(on_0_1000) proj <- round(proj*1000)
  
  return(proj)
}





# GAM Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('GAM_biomod2_model',
         representation(model_subclass='character'),
         contains = 'biomod2_model',
         prototype(model_class = 'GAM',
                   model_subclass = 'GAM_mgcv'),
         validity = function(object){
           # check model class
           if(! (object@model_subclass %in% c('GAM_mgcv', 'GAM_gam', 'BAM_mgcv') ) )
             return(FALSE)
           
           if(object@model_subclass %in% c('GAM_mgcv','GAM_gam')) 
             if(sum(! ( c("gam","glm","lm") %in% class(object@model) ) ) > 0)
               return(FALSE)
           
           if(object@model_subclass == 'BAM_mgcv') 
             if(sum(! ( c("bam","gam","glm","lm") %in% class(object@model) ) ) > 0)
               return(FALSE)           
             
           return(TRUE)
         })

setMethod('predict', signature(object = 'GAM_biomod2_model'),
          function(object, newdata, ...){
            
            args <- list(...)
            if(object@model_subclass %in% c("GAM_mgcv","BAM_mgcv")){
              if( ! ("package:mgcv" %in% search()) ){ require(mgcv,quietly=TRUE) }
            }
            
            if(object@model_subclass == "GAM_gam"){
              if( ! ("package:gam" %in% search()) ){ require(gam,quietly=TRUE) }
            }
            
            if(inherits(newdata, 'Raster')){            
              return(.predict.GAM_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.GAM_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }
            
          })

.predict.GAM_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000
  
  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  
  proj <- .testnull(object = getFormalModel(object), Prev = 0.5 , dat = newdata)
  
  if(length(getScalingModel(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = getScalingModel(object), Prev = 0.5 , dat = proj)
  }
  
  if(on_0_1000) proj <- round(proj*1000)
  
  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    writeRaster(proj, filename=filename, overwrite=overwrite)
    proj <- raster(filename)
  }
  
  return(proj)
}

.predict.GAM_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  
  proj <- as.numeric(.testnull(object = getFormalModel(object), Prev = 0.5 , dat = as.data.frame(newdata)))

  if(length(getScalingModel(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = getScalingModel(object), Prev = 0.5 , dat = proj)
  }
  
  if(on_0_1000) proj <- round(proj*1000)
  
  return(proj)
}





# GBM Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('GBM_biomod2_model',
         representation(n.trees_optim = 'numeric'),
         contains = 'biomod2_model',
         prototype(model_class = 'GBM',
                   n.trees_optim = 1000),
         validity = function(object){
           # check model class
           if(sum(! ( c("gbm") %in% class(object@model) ) ) > 0) return(FALSE)
           return(TRUE)
         })

setMethod('predict', signature(object = 'GBM_biomod2_model'),
          function(object, newdata, ...){
            
            args <- list(...)
            
            if( ! ("package:gbm" %in% search()) ){ require(gbm,quietly=TRUE) }
            
            if(inherits(newdata, 'Raster')){            
              return(.predict.GBM_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.GBM_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }
            
          })

.predict.GBM_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000
  
  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  
  proj <- predict(newdata, model=getFormalModel(object), n.trees = object@n.trees_optim, type = "response")
  
  if(length(getScalingModel(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = getScalingModel(object), Prev = 0.5 , dat = proj)
  }
  
  if(on_0_1000) proj <- round(proj*1000)
  
  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    writeRaster(proj, filename=filename, overwrite=overwrite)
    proj <- raster(filename)
  }
  
  return(proj)
}

.predict.GBM_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  
  proj <- as.numeric(predict(getFormalModel(object), as.data.frame(newdata), n.trees = object@n.trees_optim, type = "response"))
  
  if(length(getScalingModel(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = getScalingModel(object), Prev = 0.5 , dat = proj)
  }
  
  if(on_0_1000) proj <- round(proj*1000)
  
  return(proj)
}





# GLM Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('GLM_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype(model_class = 'GLM'),
         validity = function(object){
           # check model class
           if(sum(! ( c("glm", "lm") %in% class(object@model) ) ) > 0) return(FALSE)
           return(TRUE)
         })

setMethod('predict', signature(object = 'GLM_biomod2_model'),
          function(object, newdata, ...){
            
            args <- list(...)
            
            if( ! ("package:stats" %in% search()) ){ require(stats,quietly=TRUE) }
            
            if(inherits(newdata, 'Raster')){            
              return(.predict.GLM_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.GLM_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }
            
          })

.predict.GLM_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000
  
  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  
  proj <- .testnull(object = getFormalModel(object), Prev = 0.5 , dat = newdata)
  
  if(length(getScalingModel(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = getScalingModel(object), Prev = 0.5 , dat = proj)
  }
  
  if(on_0_1000) proj <- round(proj*1000)
  
  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    writeRaster(proj, filename=filename, overwrite=overwrite)
    proj <- raster(filename)
  }
  
  return(proj)
}

.predict.GLM_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  
  proj <- as.numeric(.testnull(object = getFormalModel(object), Prev = 0.5 , dat = as.data.frame(newdata)))
  
  if(length(getScalingModel(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = getScalingModel(object), Prev = 0.5 , dat = proj)
  }
  
  if(on_0_1000) proj <- round(proj*1000)
  
  return(proj)
}





# MARS Class =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('MARS_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype(model_class = 'MARS'),
         validity = function(object){
           # check model class
           if(sum(! ( c("mars") %in% class(object@model) ) ) > 0) return(FALSE)
           return(TRUE)
         })

setMethod('predict', signature(object = 'MARS_biomod2_model'),
          function(object, newdata, ...){
            
            args <- list(...)
            
            if( ! ("package:mda" %in% search()) ){ require(mda,quietly=TRUE) }
            
            if(inherits(newdata, 'Raster')){            
              return(.predict.MARS_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.MARS_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }
            
          })

.predict.MARS_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000
  
  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  
  proj <- predict(newdata, model=getFormalModel(object))
  
  if(length(getScalingModel(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = getScalingModel(object), Prev = 0.5 , dat = proj)
  }
  
  if(on_0_1000) proj <- round(proj*1000)
  
  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    writeRaster(proj, filename=filename, overwrite=overwrite)
    proj <- raster(filename)
  }
  
  return(proj)
}

.predict.MARS_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  
  proj <- as.numeric(predict(getFormalModel(object), as.data.frame(newdata)))
  
  if(length(getScalingModel(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = getScalingModel(object), Prev = 0.5 , dat = proj)
  }
  
  if(on_0_1000) proj <- round(proj*1000)
  
  return(proj)
}





# MAXENT Class =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('MAXENT_biomod2_model',
         representation(model_output_dir = 'character'),
         contains = 'biomod2_model',
         prototype(model_class = 'MAXENT'),
         validity = function(object){
           # check model class
#            if(sum(! ( c("randomForest.formula", "randomForest") %in% class(object@model) ) ) > 0) return(FALSE)
           return(TRUE)
         })

setMethod('predict', signature(object = 'MAXENT_biomod2_model'),
          function(object, newdata, ...){
            
            args <- list(...)
            
            if(inherits(newdata, 'Raster')){            
              return(.predict.MAXENT_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.MAXENT_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }
            
          })


.predict.MAXENT_biomod2_model.RasterStack <- function(object, newdata,  ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000
  rm_tmp_files <- args$rm_tmp_files
  temp_workdir <- args$temp_workdir
  
  if (is.null(temp_workdir)) temp_workdir <- paste("maxentWDtmp", format(Sys.time(), "%s"), sep="")
  if (is.null(rm_tmp_files)) rm_tmp_files <- TRUE
  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  
  .Prepare.Maxent.Proj.WorkDir(Data = newdata, proj_name = file.path(object@resp_name,temp_workdir))
  
  cat("\n\t\tRunning Maxent...")
  system(command=paste("java -cp ", file.path(object@model_options$path_to_maxent.jar, "maxent.jar"),
                       " density.Project \"", 
                       file.path(object@model_output_dir, sub("_MAXENT",".lambdas",object@model_name, fixed=T)),"\" ", 
                       file.path(object@resp_name, temp_workdir, "MaxentTmpData","Proj"), " ", 
                       file.path(object@resp_name, temp_workdir, "MaxentTmpData", "projMaxent.grd") , 
                       " doclamp=false visible=false autorun nowarnings notooltips", sep=""), wait = TRUE)
  
  cat("\n\t\tReading Maxent outputs...")
  proj <- raster(file.path(object@resp_name, temp_workdir , "MaxentTmpData","projMaxent.grd"))
  
  if(length(getScalingModel(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = getScalingModel(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)
  
  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    writeRaster(proj, filename=filename, overwrite=overwrite)
    proj <- raster(filename)
  } else if(!inMemory(proj)){
    proj <- readAll(proj) # to prevent from tmp files removing
  }
  
  if(!is.null(rm_tmp_files)){
    if(rm_tmp_files){
      unlink(x=file.path(object@resp_name, temp_workdir), recursive=TRUE, force=TRUE )
    }
  }
  
  return(proj)
}

.predict.MAXENT_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  temp_workdir <- args$temp_workdir
  rm_tmp_files <- args$rm_tmp_files
  xy <- args$xy
  
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  if (is.null(temp_workdir)) temp_workdir <- paste("maxentWDtmp", format(Sys.time(), "%s"), sep="")
  if (is.null(rm_tmp_files)) rm_tmp_files <- TRUE
  
#   if( is.null(xy) ){
#     if( sum(c('x','y') %in% colnames(newdata) ) == 2 ){
#       coor_col <- c( which(colnames(newdata) == 'x'), which(colnames(newdata) == 'y') )
#       xy <- newdata[,coor_col]
#       newdata <- newdata[,- coor_col]
#     } else { 
#       xy <- data.frame(x=rep(0,nrow(newdata)), y=rep(0,nrow(newdata)))
#     }
#   }
  
  ## no xy needed for models projections
  xy <- NULL
  
  .Prepare.Maxent.Proj.WorkDir(Data = as.data.frame(newdata), xy = xy , proj_name = file.path(object@resp_name,temp_workdir))
  
  cat("\n\t\tRunning Maxent...")
  system(command=paste("java -cp ", file.path(object@model_options$path_to_maxent.jar, "maxent.jar"),
                       " density.Project \"", 
                       file.path(object@model_output_dir, sub("_MAXENT",".lambdas",object@model_name, fixed=T)),"\" ", 
                       file.path(object@resp_name, temp_workdir, "MaxentTmpData","Proj_swd.csv"), " ", 
                       file.path(object@resp_name, temp_workdir, "MaxentTmpData", "projMaxent.asc") , 
                       " doclamp=false", sep=""), wait = TRUE)
  
  cat("\n\t\tReading Maxent outputs...")
  proj <- as.numeric(read.asciigrid(file.path(object@resp_name, temp_workdir , "MaxentTmpData", "projMaxent.asc"))@data[,1])
  
  if(!is.null(rm_tmp_files)){
    if(rm_tmp_files){
      unlink(file.path(object@resp_name, temp_workdir),recursive=TRUE,force=TRUE)
    }
  }

  if(length(getScalingModel(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = getScalingModel(object), Prev = 0.5 , dat = proj)
  }
  
  if(on_0_1000) proj <- round(proj*1000)
  
  return(proj)
  
}





# RF Class =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('RF_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype(model_class = 'RF'),
         validity = function(object){
           # check model class
           if(sum(! ( c("randomForest.formula", "randomForest") %in% class(object@model) ) ) > 0) return(FALSE)
           return(TRUE)
         })

setMethod('predict', signature(object = 'RF_biomod2_model'),
          function(object, newdata, ...){
            
            args <- list(...)
            
            if( ! ("package:randomForest" %in% search()) ){ require(randomForest,quietly=TRUE) }
            
            if(inherits(newdata, 'Raster')){            
              return(.predict.RF_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.RF_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }
            
          })

.predict.RF_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000
  
  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  
  proj <- predict(newdata, model=getFormalModel(object), type='prob', index=2)
  
  if(length(getScalingModel(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = getScalingModel(object), Prev = 0.5 , dat = proj)
  }
  
  if(on_0_1000) proj <- round(proj*1000)
  
  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    writeRaster(proj, filename=filename, overwrite=overwrite)
    proj <- raster(filename)
  }
  
  return(proj)
}

.predict.RF_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  
  proj <- as.numeric(predict(getFormalModel(object), as.data.frame(newdata), type='prob')[,'1'])

  if(length(getScalingModel(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = getScalingModel(object), Prev = 0.5 , dat = proj)
  }
  
  if(on_0_1000) proj <- round(proj*1000)
  
  return(proj)
}





# SRE Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('SRE_biomod2_model',
         representation(extremal_conditions='data.frame'),
         contains = 'biomod2_model',
         prototype(model_class = 'SRE'),
         validity = function(object){
           return(TRUE)
         })

setMethod('predict', signature(object = 'SRE_biomod2_model'),
          function(object, newdata, ...){
            
            args <- list(...)
            
            if(inherits(newdata, 'Raster')){            
              return(.predict.SRE_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.SRE_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }
            
          })

.predict.SRE_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000
  
  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  
  proj <- .sre.projection(NewData=newdata, ExtremCond=object@extremal_conditions)
  
  if(on_0_1000) proj <- round(proj*1000)

  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    writeRaster(proj, filename=filename, overwrite=overwrite)
    proj <- raster(filename)
  }
  
  return(proj)
}

.predict.SRE_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  
  proj <- .sre.projection(NewData=newdata, ExtremCond=object@extremal_conditions)
  
  if(on_0_1000) proj <- round(proj*1000)
  
  return(proj)
}

