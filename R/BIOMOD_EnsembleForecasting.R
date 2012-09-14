'BIOMOD_EnsembleForecasting' <- function( projection.output,
                                          EM.output,
                                          total.consensus = FALSE,
                                          binary.meth = NULL,
                                          filtered.meth = NULL ){
  .bmCat("Do Ensemble Models Projections")
  
#   cat("\n*** Dennis Debugged Version ***\n")
  # 1. args checking -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  args <- .BIOMOD_EnsembleForecasting.check.args( projection.output,
                                                  EM.output,
                                                  total.consensus,
                                                  binary.meth,
                                                  filtered.meth )
  
  total.consensus <- args$total.consensus
  rm('args')
                                                  
  
  # 2. Do the ensemble modeling
  for( em.comp in EM.output@em.computed){
    cat("\n\n***", em.comp, "...")
    cleanWD <- c(ls(),"cleanWD")
    for( em.algo in getEMalgos(EM.output, em.comp) ){
      cat("\n   >", em.algo)
      
      # 1. Mean of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
      if (em.algo == 'em.mean'){
        if(projection.output@type == 'RasterStack'){
          ef.mean <- raster:::mean(raster:::subset(getProjection(projection.output), 
                                                 getEMkeptModels(EM.output, em.comp)))
        } else if(projection.output@type == 'array'){
          ef.mean <- rowMeans(getProjection(projection.output, as.data.frame = TRUE)[,getEMkeptModels(EM.output, em.comp)])
        } else if(projection.output@type == 'character'){ 
          # get models to load
          projToLoad <- c()
          for(mod in getEMkeptModels(EM.output, em.comp)){
            projToLoad <- c(projToLoad, grep(mod, projection.output@proj@val, fixed=T, value=TRUE))
          }
          if(length(projToLoad)<1){
            cat("\nnot done because of invalid models names!")
          } else{
            # load the first raster (as mask)
            ef.mean <- get(load(paste(projection.output@proj@link, projToLoad[1], sep="")))
            rm(list=paste(projToLoad[1]))
            # sum all projections
            if(length(projToLoad) > 1 ){
              for(ptl in projToLoad[-1]){
                ef.mean <- ef.mean + get(load(paste(projection.output@proj@link, ptl, sep="")))
                rm(list=paste(ptl))
              }
            }
            # dividing by number of projection to get mean
            ef.mean <- ef.mean / length(projToLoad)
          }
        } else{
          cat("Unsupported yet !")
        }
      }
      
      # 2. CV of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
      if(em.algo == 'em.cv'){
        if(projection.output@type == 'RasterStack'){
          ef.cv <- round(raster:::cv(raster:::subset(getProjection(projection.output), 
                                                 getEMkeptModels(EM.output, em.comp))), digits=2)
        } else if(projection.output@type == 'array'){
          ef.sd <- apply(getProjection(projection.output, as.data.frame = TRUE)[,getEMkeptModels(EM.output, em.comp)],1,sd)
          if(!exists('ef.mean')){
            ef.mean <- rowMeans(getProjection(projection.output, as.data.frame = TRUE)[,getEMkeptModels(EM.output, em.comp)])
          }
          ef.cv <- round(ef.sd / ef.mean, 2)
          # putting to 0 points where mean = 0
          ef.cv[ ef.mean == 0 ] <- 0
        } else if(projection.output@type == 'character'){ 
          # get models to load
          projToLoad <- c()
          for(mod in getEMkeptModels(EM.output, em.comp)){
            projToLoad <- c(projToLoad, grep(mod, projection.output@proj@val, fixed=T, value=TRUE))
          }
          if(length(projToLoad)<1){
            cat("\nnot done because of invalid models names!")
          } else{
            
            if(!exists('ef.mean')){
              # load the first raster (as mask)
              ef.mean <- get(load(paste(projection.output@proj@link, projToLoad[1], sep="")))
              rm(list=paste(projToLoad[1]))
              # sum all projections
              if(length(projToLoad) > 1 ){
                for(ptl in projToLoad[-1]){
                  ef.mean <- ef.mean + get(load(paste(projection.output@proj@link, ptl, sep="")))
                  rm(list=paste(ptl))
                }
              }
              # dividing by number of projection to get mean
              ef.mean <- ef.mean / length(projToLoad)              
            }
            
            if(!exists('ef.sd')){
              # load the first raster (as mask)
              ef.sd <- (get(load(paste(projection.output@proj@link, projToLoad[1], sep=""))) - ef.mean)^2
              rm(list=paste(projToLoad[1]))
              # sum all projections
              if(length(projToLoad) > 1 ){
                for(ptl in projToLoad[-1]){
                  ef.sd <- ef.sd + (get(load(paste(projection.output@proj@link, ptl, sep=""))) - ef.mean)^2
                  rm(list=paste(ptl))
                }
              }
              # dividing by number of projection to get mean and keepin the square root to have sd
              ef.sd <- sqrt(ef.sd / length(projToLoad))
            }
            
            ef.cv <- round(ef.sd / ef.mean, 2)
            # putting to 0 points where mean = 0
            ef.cv[ ef.mean[] == 0 ] <- 0            
          }
        } else { 
          cat("Unsupported yet !")
        }        
      }
      
      # 3. Median of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
      if(em.algo == 'em.median'){
        if(projection.output@type == 'RasterStack'){
          ef.median <- round(calc(raster:::subset(getProjection(projection.output), 
                                                 getEMkeptModels(EM.output, em.comp)), median))
        } else if(projection.output@type == 'array'){
        ef.median <- round(apply(getProjection(projection.output, as.data.frame = TRUE)[,getEMkeptModels(EM.output, em.comp)], 1, median))
        } else if(projection.output@type == 'character'){
          cat("\n      ! Not done because projection RasterStack seems to be too heavy !")
        } else { 
          cat("Unsupported yet !")
        }
      }
      
      # 4. CI of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
      if(em.algo == 'em.ci.inf'){
        if(projection.output@type == 'RasterStack'){
          if(!exists('ef.mean')){
            ef.mean <- raster:::mean(raster:::subset(getProjection(projection.output), getEMkeptModels(EM.output, em.comp)))
          }
          if(!exists('ef.sd')){
            ef.sd <- calc(raster:::subset(getProjection(projection.output), getEMkeptModels(EM.output, em.comp)), sd)
          }
          ef.ci.inf <- round(ef.mean - qt(1-EM.output@em.ci.alpha/2, df = length(getEMkeptModels(EM.output, em.comp)) + 1 ) / sqrt(length(getEMkeptModels(EM.output, em.comp))) * ef.sd)
          ef.ci.inf <- reclassify(ef.ci.inf, c(-Inf,0,0))
        } else if(projection.output@type == 'array'){
          if(!exists('ef.mean')){
            ef.mean <- mean(getProjection(projection.output, as.data.frame = TRUE)[,getEMkeptModels(EM.output, em.comp)])
          }
          if(!exists('ef.sd')){
            ef.sd <- apply(getProjection(projection.output, as.data.frame = TRUE)[,getEMkeptModels(EM.output, em.comp)],1,sd)
          }
          ef.ci.inf <- round(ef.mean - qt(1-EM.output@em.ci.alpha/2, df = length(getEMkeptModels(EM.output, em.comp)) + 1 ) / sqrt(length(getEMkeptModels(EM.output, em.comp))) * ef.sd)
          ef.ci.inf[ef.ci.inf<0] <- 0
        } else if(projection.output@type == 'character'){
          # get models to load
          projToLoad <- c()
          for(mod in getEMkeptModels(EM.output, em.comp)){
            projToLoad <- c(projToLoad, grep(mod, projection.output@proj@val, fixed=T, value=TRUE))
          }
          if(length(projToLoad)<1){
            cat("\nnot done because of invalid models names!")
          } else {
            
            if(!exists('ef.mean')){
              # load the first raster (as mask)
              ef.mean <- get(load(paste(projection.output@proj@link, projToLoad[1], sep="")))
              rm(list=paste(projToLoad[1]))
              # sum all projections
              if(length(projToLoad) > 1 ){
                for(ptl in projToLoad[-1]){
                  ef.mean <- ef.mean + get(load(paste(projection.output@proj@link, ptl, sep="")))
                  rm(list=paste(ptl))
                }
              }
              # dividing by number of projection to get mean
              ef.mean <- ef.mean / length(projToLoad)              
            }
            
            if(!exists('ef.sd')){
              # load the first raster (as mask)
              ef.sd <- (get(load(paste(projection.output@proj@link, projToLoad[1], sep=""))) - ef.mean)^2
              rm(list=paste(projToLoad[1]))
              # sum all projections
              if(length(projToLoad) > 1 ){
                for(ptl in projToLoad[-1]){
                  ef.sd <- ef.sd + (get(load(paste(projection.output@proj@link, ptl, sep=""))) - ef.mean)^2
                  rm(list=paste(ptl))
                }
              }
              # dividing by number of projection to get mean and keepin the square root to have sd
              ef.sd <- sqrt(ef.sd / length(projToLoad))
            }
            
            ef.ci.inf <- round(ef.mean - qt(1-EM.output@em.ci.alpha/2, df = length(projToLoad) + 1 ) / sqrt(length(projToLoad)) * ef.sd)
            ef.ci.inf <- reclassify(ef.ci.inf, c(-Inf,0,0))
          }             
        } else { 
          cat("Unsupported yet !")
        }
      }
      
      if(em.algo == 'em.ci.sup'){
        if(projection.output@type == 'RasterStack'){
          if(!exists('ef.mean')){
            ef.mean <- raster:::mean(raster:::subset(getProjection(projection.output), getEMkeptModels(EM.output, em.comp)))
          }
          if(!exists('ef.sd')){
            ef.sd <- calc(raster:::subset(getProjection(projection.output), getEMkeptModels(EM.output, em.comp)), sd)
          }
          ef.ci.sup <- round(ef.mean + qt(1-EM.output@em.ci.alpha/2, df = length(getEMkeptModels(EM.output, em.comp)) + 1 ) / sqrt(length(getEMkeptModels(EM.output, em.comp))) * ef.sd)
          ef.ci.sup <- reclassify(ef.ci.sup, c(1000,+Inf,1000))
        } else if(projection.output@type == 'array'){
          if(!exists('ef.mean')){
            ef.mean <- mean(getProjection(projection.output, as.data.frame = TRUE)[,getEMkeptModels(EM.output, em.comp)])
          }
          if(!exists('ef.sd')){
            ef.sd <- apply(getProjection(projection.output, as.data.frame = TRUE)[,getEMkeptModels(EM.output, em.comp)],1,sd)
          }
          ef.ci.sup <- round(ef.mean + qt(1-EM.output@em.ci.alpha/2, df = length(getEMkeptModels(EM.output, em.comp)) + 1 ) / sqrt(length(getEMkeptModels(EM.output, em.comp))) * ef.sd)
          ef.ci.inf[ef.ci.inf>1000] <- 1000
        } else if(projection.output@type == 'character'){
          # get models to load
          projToLoad <- c()
          for(mod in getEMkeptModels(EM.output, em.comp)){
            projToLoad <- c(projToLoad, grep(mod, projection.output@proj@val, fixed=T, value=TRUE))
          }
          if(length(projToLoad)<1){
            cat("\nnot done because of invalid models names!")
          } else {
            
            if(!exists('ef.mean')){
              # load the first raster (as mask)
              ef.mean <- get(load(paste(projection.output@proj@link, projToLoad[1], sep="")))
              rm(list=paste(projToLoad[1]))
              # sum all projections
              if(length(projToLoad) > 1 ){
                for(ptl in projToLoad[-1]){
                  ef.mean <- ef.mean + get(load(paste(projection.output@proj@link, ptl, sep="")))
                  rm(list=paste(ptl))
                }
              }
              # dividing by number of projection to get mean
              ef.mean <- ef.mean / length(projToLoad)              
            }
            
            if(!exists('ef.sd')){
              # load the first raster (as mask)
              ef.sd <- (get(load(paste(projection.output@proj@link, projToLoad[1], sep=""))) - ef.mean)^2
              rm(list=paste(projToLoad[1]))
              # sum all projections
              if(length(projToLoad) > 1 ){
                for(ptl in projToLoad[-1]){
                  ef.sd <- ef.sd + (get(load(paste(projection.output@proj@link, ptl, sep=""))) - ef.mean)^2
                  rm(list=paste(ptl))
                }
              }
              # dividing by number of projection to get mean and keepin the square root to have sd
              ef.sd <- sqrt(ef.sd / length(projToLoad))
            }
            
            ef.ci.sup <- round(ef.mean + qt(1-EM.output@em.ci.alpha/2, df = length(projToLoad) + 1 ) / sqrt(length(projToLoad)) * ef.sd)
            ef.ci.sup <- reclassify(ef.ci.sup, c(1000,+Inf,1000))
          }             
        } else { 
          cat("Unsupported yet !")
        }
      }

      # 5. Committee averaging of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
      if(em.algo == 'em.ca'){
        ### May be a good idea to take the binary projection if computed
        models.kept.tresh <- eval(parse(text = paste("EM.output@em.bin.tresh$", em.comp, sep="")))
      
        if(projection.output@type == 'RasterStack'){
          ef.ca <- round(raster:::mean(BinaryTransformation(raster:::subset(getProjection(projection.output), 
                                                 getEMkeptModels(EM.output, em.comp)), models.kept.tresh)) * 1000)
        } else if(projection.output@type == 'array'){
          ef.ca <- round(apply(as.data.frame(BinaryTransformation(getProjection(projection.output, as.data.frame = TRUE)[,getEMkeptModels(EM.output, em.comp)] ,models.kept.tresh)), 1, mean)*1000)
        } else if(projection.output@type == 'character'){
          # get models to load
          projToLoad <- c()
          for(mod in getEMkeptModels(EM.output, em.comp)){
            projToLoad <- c(projToLoad, grep(mod, projection.output@proj@val, fixed=T, value=TRUE))
          }
          if(length(projToLoad)<1){
            cat("\nnot done because of invalid models names!")
          } else {
            # load the first raster (as mask)
            ef.ca <- BinaryTransformation(get(load(paste(projection.output@proj@link, projToLoad[1], sep=""))),  models.kept.tresh[1])
            rm(list=paste(projToLoad[1]))
            # sum all projections
            if(length(projToLoad) > 1 ){
              for(ptl in projToLoad[-1]){
                ef.ca <- ef.ca + BinaryTransformation(get(load(paste(projection.output@proj@link, ptl, sep=""))), 
                                                      models.kept.tresh[which(projToLoad == ptl)])
                rm(list=paste(ptl))
              }
            }
            # dividing by number of projection to get mean
            ef.ca <- ef.ca / length(projToLoad)            
          }
        } else { 
          cat("Unsupported yet !")
        }
      }
      
      # 6. weighted mean of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
      if(em.algo == 'em.pmw'){
#         cat("\n   > Prababilities wegthing mean...")
#         
#         # remove SRE models if ROC
#         if(eval.m == 'ROC'){
#           sre.id <- grep("_SRE", models.kept)
#           if(length(sre.id)>0){
#             cat("\n      ! SRE modeling were switched off")
#             models.kept <- models.kept[-sre.id]
#             models.kept.scores <- models.kept.scores[-sre.id]
#             prediction.kept <- prediction.kept[,models.kept]
#           }
#         }
        
        models.kept.scores = EM.output@em.weight[[em.comp]]
      	### Compute ensemble forecast
        if(projection.output@type == 'RasterStack'){
          if(length(getEMkeptModels(EM.output, em.comp)) > 1){
            ef.pmw <- round(sum(raster:::subset(getProjection(projection.output), getEMkeptModels(EM.output, em.comp)) * models.kept.scores))
          } else{
            ef.pmw <- raster:::subset(getProjection(projection.output), getEMkeptModels(EM.output, em.comp))
          }

        } else if(projection.output@type == 'array'){
          if(length(getEMkeptModels(EM.output, em.comp)) > 1){
            ef.pmw <- round(as.vector(data.matrix(getProjection(projection.output, as.data.frame = TRUE)[,getEMkeptModels(EM.output, em.comp)])%*% models.kept.scores))
          } else{
            ef.pmw <- as.vector(getProjection(projection.output, as.data.frame = TRUE)[,getEMkeptModels(EM.output, em.comp)])
          }
        } else if(projection.output@type == 'character'){
          # get models to load
          projToLoad <- c()
          for(mod in getEMkeptModels(EM.output, em.comp)){
            projToLoad <- c(projToLoad, grep(mod, projection.output@proj@val, fixed=T, value=TRUE))
          }
          if(length(projToLoad)<1){
            cat("\nnot done because of invalid models names!")
          } else {
            # load the first raster (as mask)
            ef.pmw <- get(load(paste(projection.output@proj@link, projToLoad[1], sep=""))) * models.kept.scores[1]
            rm(list=paste(projToLoad[1]))
            # sum all projections
            if(length(projToLoad) > 1 ){
              for(ptl in projToLoad[-1]){
                ef.pmw <- ef.pmw + get(load(paste(projection.output@proj@link, ptl, sep=""))) * models.kept.scores[which(projToLoad == ptl)]
                rm(list=paste(ptl))
              }
            }          
          }
        } else { 
          cat("Unsupported yet !")
        }
      }     
    }
    
    # 7. Assembling all computed models  and save projections on hard drive -=-=-=-=-=-=-=-=-=-=-= #
    ef.potential <- c('ef.mean','ef.cv','ef.ci.inf','ef.ci.sup','ef.median','ef.ca','ef.pmw')
    ef.computed <- ef.potential[unlist(lapply(ef.potential, exists, envir = environment()))]
    
    if(projection.output@type == 'RasterStack' | projection.output@type == 'character'){
      eval(parse(text = paste(em.comp ,"<- raster:::stack(", toString(ef.computed), ")", sep="")))
      eval(parse(text = paste("names(", em.comp, ") <-  ef.computed", sep="")))
    } else if(projection.output@type == 'array'){
      eval(parse(text = paste(em.comp ," <- cbind(", toString(ef.computed), ")", sep="")))
    } else{
      cat("Unsupported yet !")
    }
    
    eval(parse(text = paste("save(", em.comp, ", file = '", 
                            projection.output@sp.name, .Platform$file.sep, "proj_",
                            projection.output@proj.names, .Platform$file.sep,
                            em.comp, "')", sep="")))
    
    # 8. Doing Binary Transformation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    if(!is.null(binary.meth)){
      for(bin.meth in binary.meth){
        cuts <- getEMeval(EM.output, em.comp)[[1]][bin.meth, "Cutoff", ]
        cuts <- cuts[sub("ef.","em.", ef.computed)]
#         cat("\n*** cuts selected !")
#         cat("\n*** cuts =", cuts)
        if(sum(is.na(cuts)) > 0) {
          warning(em.comp, " : ",toString(ef.computed[which(is.na(cuts))])," cuts was automaticly set to 500.. That can lead to strange binary results.") 
          cuts[which(is.na(cuts))] <- 500
          }

        eval(parse(text=paste(em.comp,".bin.", bin.meth," <- BinaryTransformation(", em.comp,", cuts)" , sep="")))
          
        if(projection.output@type == 'RasterStack' | projection.output@type == 'character'){
          eval(parse(text=paste("names(", em.comp,".bin.", bin.meth ,") <- paste(names(",em.comp,"), '.bin',sep='')", sep="")))
        }
        eval(parse(text = paste("save(", em.comp,".bin.", bin.meth, ", file = '", 
                                projection.output@sp.name, .Platform$file.sep, "proj_",
                                projection.output@proj.names, .Platform$file.sep,
                                em.comp,".bin.", bin.meth, "')", sep="")))
      }
    }
    
    # 9. Doing Filtred Transformation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    if(!is.null(filtered.meth)){
      cat("\n      ! Filtered transformation not supported yet!")
    }

    # 10. Cleaning work dir -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    rm(list = ls()[!(ls() %in% cleanWD)])
  }

  # Doing Total Consensus Projection -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  ## NOTE : We have to think about the best way to make this concencus.. At time, let's simply do a mean of EF conputed
  if(total.consensus){
    cat("\n\n*** Projections Consensus...")
#     if(projection.output@type == 'RasterStack'){
#       ef.cons <- raster:::stack()
#     } else if(projection.output@type == 'array'){
#       ef.cons <- array()
#     } else { 
#       cat("Unsupported yet !")
#     }
    ### Initialisation
    ef.cons <- get(load(paste(projection.output@sp.name, .Platform$file.sep, "proj_", projection.output@proj.names,
                   .Platform$file.sep, EM.output@em.computed[1], sep="")))
    
    if(class(ef.cons) == 'RasterStack'){
      ef.computed <- names(ef.cons)
    } else{
      ef.computed <- colnames(ef.cons)
    }

    rm(list = paste(EM.output@em.computed[1], sep=""))
    
    ### Filling
    for(em.comp in EM.output@em.computed[-1]){
      ef.tmp <- get(load(paste(projection.output@sp.name, .Platform$file.sep, "proj_", projection.output@proj.names,
                 .Platform$file.sep, em.comp, sep="")))
      rm(list = paste(em.comp, sep=""))
      ef.cons <- ef.cons + ef.tmp
    }
    
    eval(parse(text=paste(projection.output@sp.name, "_TotalConsensus <- ef.cons / length(EM.output@em.computed)",sep=""))) 
    rm(list = c('ef.tmp', 'ef.cons'))
      
    ### Saving
    eval(parse(text = paste("save(", projection.output@sp.name, "_TotalConsensus, file = '", 
                            projection.output@sp.name, .Platform$file.sep, "proj_",
                            projection.output@proj.names, .Platform$file.sep,
                            projection.output@sp.name, "_TotalConsensus')", sep="")))
    ### Bin Trans
    if(!is.null(binary.meth)){
      for(bin.meth in binary.meth){
        cuts <- t(as.data.frame(sapply(EM.output@em.computed, function(em.comp){
          return(getEMeval(EM.output, em.comp)[[1]][bin.meth, "Cutoff", ])})))
        cuts <- apply(cuts,2,mean, na.rm=T)
        cuts <- cuts[sub("ef.", "em.", ef.computed)]
    
        eval(parse(text=paste(projection.output@sp.name, "_TotalConsensus.bin.", bin.meth," <- BinaryTransformation(",
                              projection.output@sp.name, "_TotalConsensus, cuts)" , sep="")))
        if(projection.output@type == 'RasterStack' | projection.output@type == 'character'){
          eval(parse(text=paste("names(", projection.output@sp.name, "_TotalConsensus.bin.", bin.meth,") <- paste(names(",projection.output@sp.name, "_TotalConsensus), '.bin',sep='')", sep="")))
        }
        eval(parse(text = paste("save(", projection.output@sp.name, "_TotalConsensus.bin.", bin.meth, ", file = '", 
                                projection.output@sp.name, .Platform$file.sep, "proj_",
                                projection.output@proj.names, .Platform$file.sep,
                                projection.output@sp.name, "_TotalConsensus.bin.", bin.meth, "')", sep="")))
      }

    }

  }
  .bmCat("Done") 
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.BIOMOD_EnsembleForecasting.check.args <- function( projection.output,
                                                    EM.output,
                                                    total.consensus,
                                                    binary.meth,
                                                    filtered.meth ){

                              
  if(total.consensus){
    if(length(EM.output@em.computed) < 2){
      cat("\n      ! Total consensus projection was switched off because only one Ensemble modeling was done")
      total.consensus <- FALSE
    }
  }
    
  if(!is.null(binary.meth)){
    if(sum(!(binary.meth %in% EM.output@eval.metric))){
      stop(paste("binary methods must be compatible with Ensemble Modeling evaluation metrics (e.g. ",
                 toString(EM.output@eval.metric)," )", sep=""))
    }
  }
  
  if(!is.null(filtered.meth)){
    if(sum(!(filtered.meth %in% EM.output@eval.metric))){
      stop(paste("filtering methods must be compatible with Ensemble Modeling evaluation metrics (e.g. ",
                 toString(EM.output@eval.metric)," )", sep=""))
    }
  }
                                                    
  return(list(projection.output = projection.output,
              EM.output = EM.output,
              total.consensus = total.consensus,
              binary.meth = binary.meth,
              filtered.meth = filtered.meth))
  
}