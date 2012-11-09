# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# Maxent models class definition
# Damien Georges - 06/11/2012
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

# Description =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
#
#   In this file we define a biomod2 specific class to consider maxent
# modelling runs as R models.
# All information required to make projections will be stored in
# this object 
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

# 1.1 Class Definition
setClass('maxent_model',
         representation(sp_name = 'character',
                        model_name = 'character',
                        model_output_dir = 'character',
                        model_options = 'list'),
         validity = function(object){ return(TRUE) } )

# 1.2 Constructors
# setGeneric('maxent_model', 
#             def = function(sp_name, ...){
#               standardGeneric( 'maxent_model' )
#             } )

maxent_model <- function(sp_name , model_name=NULL, model_output_dir=NULL, model_options=NULL){
  maxent_mod <- new('maxent_model', 
                    sp_name = sp_name,
                    model_name = model_name,
                    model_output_dir = model_output_dir,
                    model_options = model_options)
  return(maxent_mod)
}

setMethod('show', signature('maxent_model'),
          function(object){
            .bmCat("'maxent_model'")
            cat("\nsp.name = ", object@sp_name,fill=.Options$width)
            cat("\nmodel_name = ", object@model_name,fill=.Options$width)
            cat("\nmodel_output_dir = ", object@model_output_dir,fill=.Options$width)
            cat("\n")
            cat("\nmodel_options : ")
            print(object@model_options)
            .bmCat()
          })


setMethod('predict', signature(object = 'maxent_model'),
          function(object, newdata, proj_name, ...){
#             callGeneric(object)
            args <- list(...)
            
            xy <- args$xy # may be NULL
            
            if(inherits(newdata, 'Raster')){            
              return(predict.maxent_model.RasterStack(object, newdata, proj_name, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(predict.maxent_model.data.frame(object, newdata, proj_name, ... ))
            } else{ stop("invalid newdata input") }
            
          })


predict.maxent_model.data.frame <- function(object, newdata, proj_name, ...){
  args <- list(...)
  rm_tmp_files <- args$rm_tmp_files
  xy <- args$xy
  
  if( is.null(xy) ){
    if( sum(c('x','y') %in% colnames(newdata) ) == 2 ){
      coor_col <- c( which(colnames(newdata) == 'x'), which(colnames(newdata) == 'y') )
      xy <- newdata[,coor_col]
      newdata <- newdata[,- coor_col]
    } else { 
      xy <- data.frame(x=rep(0,nrow(newdata)), y=rep(0,nrow(newdata)))
    }
  }
  
  .Prepare.Maxent.Proj.WorkDir(Data = newdata, xy = xy , proj_name = file.path(object@sp_name,proj_name))
  
  cat("\n\t\tRuning Maxent...")
  system(command=paste("java -cp ", file.path(object@model_options$path_to_maxent.jar, "maxent.jar"),
                       " density.Project \"", 
                       file.path(object@model_output_dir, sub("_MAXENT",".lambdas",object@model_name, fixed=T)),"\" ", 
                       file.path(object@sp_name, proj_name, "MaxentTmpData","Proj_swd.csv"), " ", 
                       file.path(object@sp_name, proj_name, "MaxentTmpData", "projMaxent.asc") , 
                       " doclamp=false", sep=""), wait = TRUE)
  
  cat("\n\t\tReading Maxent outputs...")
  maxent.proj <- read.asciigrid(file.path(object@sp_name, proj_name , "MaxentTmpData", "projMaxent.asc"))@data

  if(!is.null(rm_tmp_files)){
    if(rm_tmp_files){
      .Delete.Maxent.WorkDir(species.name=file.path(object@sp_name, proj_name) )
    }
  }
  
  return(as.numeric(maxent.proj[,1]))
  
}


predict.maxent_model.RasterStack <- function(object, newdata, proj_name=NULL, ...){
  args <- list(...)
  rm_tmp_files <- args$rm_tmp_files
  
  .Prepare.Maxent.Proj.WorkDir(Data = newdata, proj_name = file.path(object@sp_name,proj_name))
  
  cat("\n\t\tRuning Maxent...")
  system(command=paste("java -cp ", file.path(object@model_options$path_to_maxent.jar, "maxent.jar"),
                       " density.Project \"", 
                       file.path(object@model_output_dir, sub("_MAXENT",".lambdas",object@model_name, fixed=T)),"\" ", 
                       file.path(object@sp_name, proj_name, "MaxentTmpData","Proj"), " ", 
                       file.path(object@sp_name, proj_name, "MaxentTmpData", "projMaxent.grd") , 
                       " doclamp=false", sep=""), wait = TRUE)
  
  cat("\n\t\tReading Maxent outputs...")
  maxent.proj <- raster(file.path(object@sp_name, proj_name , "MaxentTmpData","projMaxent.grd"))
  
  if(!is.null(rm_tmp_files)){
    if(rm_tmp_files){
      .Delete.Maxent.WorkDir(species.name=file.path(object@sp_name, proj_name) )
    }
  }

  return(maxent.proj)
  
}
