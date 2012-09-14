####################################################################################################
# BIOMOD_Projection
# Damien.G
# feb 2012
####################################################################################################

# AIM :
#   Project models from BIOMOD_Modeling with different explanatory variables

# INPUT :


# OUTPUT : 



# NOTE :
#   It would be nice to add done projection to input Biomod.models.object
#   .BIOMOD_Projection.check.args <- may be reorder variables if necessary

####################################################################################################
'BIOMOD_Projection' <- function(modeling.output,
                                new.env,
                                proj.name,
                                xy.new.env = NULL,
                                selected.models = 'all',
                                binary.meth = NULL,
                                filtered.meth = NULL,
                                compress = 'xz',
                                clamping.mask = TRUE,
                                ...){
  # 0. get additional args
  add.args <- list(...)
  if(!is.null(add.args$do.stack)){
    do.stack <- add.args$do.stack
  } else{
    do.stack <- TRUE
  }
  
  if(!is.null(add.args$silent)){
    silent <- add.args$silent
  } else{
    silent <- FALSE
  }
  
  if(!silent){
    .bmCat("Do Models Projections")
  }
  
  # 1. Some Inputs args checking
  args <- .BIOMOD_Projection.check.args(modeling.output,
                                        new.env,
                                        proj.name,
                                        xy.new.env,
                                        selected.models,
                                        binary.meth,
                                        filtered.meth,
                                        compress,
                                        do.stack)#, clamping.level)
  
  proj.name <- args$proj.name
  selected.models <- args$selected.models
  binary.meth <- args$binary.meth
  filtered.meth <- args$filtered.meth
  compress <- args$compress
  do.stack <- args$do.stack
  xy.new.env <- args$xy.new.env
#   clamping.level <- args$clamping.level
  
  rm(args)
  
  # 1b. Creating the outpput object
  proj_out <- new('BIOMOD.projection.out',
              proj.names = proj.name,
              sp.name =  modeling.output@sp.name,
              expl.var.names = modeling.output@expl.var.names,
              models.projected = selected.models,
              rescaled.models = modeling.output@rescal.all.models,
              xy.coord = xy.new.env)
  
  # adapting the proj slot to projection data type (e.g. rasterStack, array)
  if(!do.stack){
    proj_out@proj <- new('BIOMOD.stored.files')
  } else if(inherits(new.env, 'Raster')){
    proj_out@proj <- new('BIOMOD.stored.raster.stack')
  } else{
    proj_out@proj <- new('BIOMOD.stored.array')
  }
    
  # 1.c creating output directory
  dir.create(paste(modeling.output@sp.name, "/proj_", proj.name, sep=""), showWarnings = FALSE, recursive = TRUE, mode = "0777")
  
  # 1.c Define the clamping mask
  if(clamping.mask){
    if(!silent) cat("\n   > defining clamping mask")
    MinMax <- getModelsInputData(modeling.output,'MinMax')
    eval(parse(text = paste("clamping.mask.",proj.name," <- .build.clamping.mask(new.env, MinMax)", sep="") )) 
    eval(parse(text = paste("save(clamping.mask.",proj.name,", file = '", modeling.output@sp.name, "/proj_", proj.name, 
                              "/", proj.name,"_clamping.mask')" ,sep="")))
  }
  
  
  
  # 2. Doing projections
  proj_out@proj@val <- Projection(models.name = selected.models,
                                     modeling.work.dir = getwd(),
                                     new.env.data  = new.env,
                                     xy = xy.new.env,
                                     proj.name = proj.name,
                                     binary.proj = binary.meth,
                                     filtred.proj = filtered.meth,
                                     models.evaluation = getModelsEvaluations(modeling.output),
                                     models.options = getModelsOptions(modeling.output),
                                     compress = compress,
                                     rescaled.models = modeling.output@rescal.all.models,
                                     do.stack = do.stack)
                                       
#   # 2.b a posteriori clamping...
#   #### TO DO : make an a priori claming
#   if(!is.null(clamping.level)){
#     cat("\n   > clamping projections...")
#     if(inherits(proj_out@proj@val, 'Raster')){
#       proj_out@proj@val <- raster:::stack(reclassify(proj_out@proj@val * reclassify( (- 1 * (clampMask)), c(-0.5,0.5,1)), c(-Inf,0,clamped.value) ))
#     } else{
#       proj_out@proj@val[which(clampMask > 0),,,] <- clamped.value
#     }
#     
#   }
  
  proj_out@type <- class(proj_out@proj@val)
  if(do.stack){
    proj_out@proj@inMemory <- TRUE
    proj_out@proj@link <- paste(modeling.output@sp.name, "/proj_", proj.name, 
                                "/", proj.name, "_", modeling.output@sp.name,
                                "_", proj_out@type, sep="")    
  } else{
    proj_out@proj@inMemory <- FALSE
    proj_out@proj@link <- paste(modeling.output@sp.name, "/proj_", proj.name, 
                                "/", sep="")    
  }

  
  
  # 3. Removing Maxent Tmp Data
  if(file.exists(paste(modeling.output@sp.name, "/proj_", proj.name,'/MaxentTmpData/',sep=''))){
    .Delete.Maxent.WorkDir(paste(modeling.output@sp.name, "/proj_", proj.name, sep=""))
  }
  
  
  if(!silent) .bmCat("Done")
  # 4. Returning output
  return(proj_out)
}

####################################################################################################
### Utilities Fuctions #############################################################################
####################################################################################################

.BIOMOD_Projection.check.args <- function(modeling.output, new.env, proj.name, xy.new.env, 
                                          selected.models, binary.meth, filtered.meth,
                                          compress, do.stack){#, clamping.level){
  ## modeling.output
  if( class(modeling.output) != 'BIOMOD.models.out'){
    stop("'modeling.output' must be the result of BIOMOD_Modeling() computation")
  }
  
  ## new.env
  # NOTE : may be reorder variables if necessary
  if( !(class(new.env) %in% c('matrix', 'data.frame', 'RasterStack') )){
    stop("'new.env' must be a matrix, a data.frame or a RasterStack")
  }
  if( class(new.env) == 'RasterStack' ){
    if(sum(!(names(new.env) %in% modeling.output@expl.var.names)) > 0 ){
      stop("'new.env' layer names don't match with explanatory variables used for buiding models")
    }
  } else{
    if(sum(!(colnames(new.env) %in% modeling.output@expl.var.names)) > 0 ){
      stop("'new.env' colnames don't match with explanatory variables used for buiding models")
    }
  }
  
  ## proj.name
  # The projection Name
  if(is.null(proj.name)){
    stop("\nYou must define a name for Projection Outpus")
  } else{
    dir.create(paste(modeling.output@sp.name,'/proj_',proj.name,'/',sep=''),
               showWarnings=FALSE)
  }
  
  ## xy.new.env
  if(!is.null(xy.new.env)  & !inherits(new.env,'Raster')){
    xy.new.env = data.matrix(xy.new.env)
    if(ncol(xy.new.env) != 2 | nrow(xy.new.env) != nrow(new.env)) stop("invalid xy coordinates argument given -- dimentions mis-match !")
  } else {
    xy.new.env = matrix()
  }
  
  ## selected.models
  if(selected.models[1] == 'all'){
    selected.models <- modeling.output@models.computed
  } else{
    selected.models <- intersect(selected.models, modeling.output@models.computed)
  }
  if(length(selected.models) < 1){
    stop('No models selected')
  }
  
  # check that given models exits
  files.check <- paste(modeling.output@sp.name,'/models/',selected.models,sep='')
  not.checked.files <- c(grep('MAXENT', files.check), grep('SRE', files.check))
  if(length(not.checked.files) > 0){files.check <- files.check[-not.checked.files]}
  missing.files <- files.check[!file.exists(files.check)]
  if( length(missing.files) > 0 ){
    stop(paste("Projection files missing : ", toString(missing.files), sep=''))
    if(length(missing.files) == length(files.check)){
      stop("Impossible to find any models, migth be a problem of working directory")
    }
  }
      
  # The binaries  and filtering transformations
  if(!is.null(binary.meth) | !is.null(filtered.meth)){
    models.evaluation <- getModelsEvaluations(modeling.output)
    if(is.null(models.evaluation)){
      warning("Binary and/or Filtred transformations of projection not ran because of models
              evaluation informations missing")
    } else{
      available.evaluation <- unique(unlist(dimnames(models.evaluation)[1]))
      if(!is.null(binary.meth)){
        if(sum(!(binary.meth %in% available.evaluation)) > 0){
          warning(paste(toString(binary.meth[!(binary.meth %in% available.evaluation)]),
                        " Binary Transformation were switched off because no correspunding",
                        " evaluation method found ", sep=""))
          binary.meth <- binary.meth[binary.meth %in% available.evaluation]
        }
      }
      
      if(!is.null(filtered.meth)){
        if(sum(!(filtered.meth %in% available.evaluation)) > 0){
          warning(paste(toString(filtered.meth[!(filtered.meth %in% available.evaluation)]),
                        " Filtred Transformation were switched off because no correspunding",
                        " evaluation method found ", sep=""))          
          filtered.meth <- filtered.meth[filtered.meth %in% available.evaluation]
        }
      }
    }
  }
  
  ## compress
  if(compress == 'xz'){
    compress <- ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz')
  }
      
  ## do.stack
  if(class(new.env) != 'RasterStack'){
    do.stack <- TRUE
  } else{
    if(do.stack){
      # test if there is memory enough to work with RasterStack
      do.stack = canProcessInMemory( raster:::subset(new.env,1), 2*length(selected.models) + nlayers(new.env) )
      if (!do.stack){ 
        cat("\n   ! Results will be saved as individual RasterLayers because of a lack of memory !")
      }
    }
  }

    
#   ## clamping checking
#   if(!is.null(clamping.level)){
#     if(!is.numeric(clamping.level)){
#       stop("clamping.level must be NULL or numeric")
#     }
#     
#     # limit clamping level
#     if( clamping.level > length(modeling.output@expl.var.names)){
#       cat("\n   ! clamping.level was down to", length(modeling.output@expl.var.names))
#       clamping.level <- length(modeling.output@expl.var.names)
#     }
#     
#     if( clamping.level < 1){
#       cat("\n   ! clamping was swich off")
#       clamping.level <- NULL
#     }  
#   }
  
  
  return(list(#modeling.output = modeling.output,
              #new.env = new.env,
              proj.name = proj.name,
              xy.new.env = xy.new.env,
              selected.models = selected.models,
              binary.meth = binary.meth,
              filtered.meth = filtered.meth,
              compress = compress,
              do.stack = do.stack))#, clamping.level = clamping.level))

}

####################################################################################################
####################################################################################################

.build.clamping.mask <- function(env, MinMax){
  if(inherits(env,'Raster')){ # raster case 
    env <- raster:::stack(env)
    env <- raster:::subset(env,names(MinMax))
    
    # create an empty mask
#     clamp.mask <- reclassify( raster:::subset(env,1, drop=TRUE), c(-Inf,Inf,0) )
    clamp.mask <- raster:::subset(env,1, drop=TRUE)
    clamp.mask[!is.na(clamp.mask[])] <- 0
    
    for(e.v in names(MinMax)){
      if(!is.null(MinMax[[e.v]]$min)){ # numeric variable
        clamp.mask <- clamp.mask + ( BinaryTransformation(raster:::subset(env, e.v, drop=TRUE), MinMax[[e.v]]$max ) + 
          (1 - BinaryTransformation(raster:::subset(env, e.v, drop=TRUE), MinMax[[e.v]]$min )) )
        
      } else if(!is.null(MinMax[[e.v]]$levels)){ # factorial variable
        clamp.mask <- clamp.mask + (raster:::subset(env, e.v, drop=TRUE) %in% MinMax[[e.v]]$levels)
      }
    }
  } else if(is.data.frame(env) | is.matrix(env) | is.numeric(env)){ # matrix and data.frame case
    env <- as.data.frame(env)
    
    # create an empty mask
    clamp.mask <- rep(0,nrow(env))
    
    for(e.v in names(MinMax)){
      if(!is.null(MinMax[[e.v]]$min)){ # numeric variable
        clamp.mask <- clamp.mask + ( BinaryTransformation(env[,e.v], MinMax[[e.v]]$max ) + 
          (1 - BinaryTransformation(env[,e.v], MinMax[[e.v]]$min )) )
      } else if(!is.null(MinMax[[e.v]]$levels)){ # factorial variable
        clamp.mask <- clamp.mask + (env[,e.v] %in% MinMax[[e.v]]$levels)
      }
    }    
    
  } else{
    stop("Unsuported env arg")
  }
  
  return(clamp.mask)
  
}