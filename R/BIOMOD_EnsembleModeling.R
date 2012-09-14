'BIOMOD_EnsembleModeling' <- function( modeling.output,
                                       chosen.models = 'all',
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = 'all',
                                       eval.metric.quality.threshold = NULL,
                                       prob.mean = TRUE,
                                       prob.cv = TRUE,
                                       prob.ci = TRUE,
                                       prob.ci.alpha = 0.05,
                                       prob.median = TRUE,
                                       committee.averaging = TRUE,
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional'){
  .bmCat("Build Ensemble Models")
  # 1. args checking -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  args <- .BIOMOD_EnsembleModeling.check.args( modeling.output,
                                               chosen.models,
                                               eval.metric,
                                               eval.metric.quality.threshold,
                                               prob.mean,
                                               prob.cv,
                                               prob.ci,
                                               prob.ci.alpha,
                                               prob.median,
                                               committee.averaging,
                                               prob.mean.weight,
                                               prob.mean.weight.decay,
                                               em.by)
  
  modeling.output <- args$modeling.output
  chosen.models <- args$chosen.models
  eval.metric <- args$eval.metric
  eval.metric.quality.threshold <- args$eval.metric.quality.threshold
  prob.mean <- args$prob.mean
  prob.cv <- args$prob.cv
  prob.ci <- args$prob.ci
  prob.ci.alpha <- args$prob.ci.alpha
  prob.median <- args$prob.median
  committee.averaging <- args$committee.averaging
  prob.mean.weight <- args$prob.mean.weight
  prob.mean.weight.decay  <- args$prob.mean.weight.decay
  em.by <- args$em.by
  
  rm('args')
  
  # 1b. creating output object and begin to fill it
#   EM <- list()
  EM <- new('BIOMOD.EnsembleModeling.out',
            sp.name = modeling.output@sp.name,
            expl.var.names = modeling.output@expl.var.names,
#             models.out.obj = new('BIOMOD.stored.models.out',
#                                  inMemory = FALSE,
#                                  link = paste(modeling.output@sp.name,"/",modeling.output@sp.name,".models.out",sep="")),
            eval.metric = eval.metric,
            eval.metric.quality.threshold = eval.metric.quality.threshold,
            em.ci.alpha = prob.ci.alpha)
  
  EM@models.out.obj@link <- paste(modeling.output@sp.name,"/",modeling.output@sp.name,".models.out",sep="")
  
  # 2. doing Ensemble modeling
  
  ## 2.1 make a list of models names that will be combined together according to by argument.
  
  # =-=-=-=-=-=-=-=- em.models.assembling function -=-=-=-=-=-=-=- #
  em.models.assembling <- function(chosen.models, em.by){
    assembl.list = list()
    
    if(em.by == 'PA_dataset'){
      for(dat in .extractModelNamesInfo(chosen.models, info='data.set')){
        assembl.list[[paste(dat,"_AllRun", sep="")]] <- chosen.models[grep(paste("_",dat,"_",sep=""), chosen.models)]
      }
      return(assembl.list)
    }
    
    if(em.by == 'algo'){
      for(algo in .extractModelNamesInfo(chosen.models, info='models')){
        assembl.list[[paste(algo,"_AllRun", sep="")]] <- chosen.models[grep(paste("_",algo,sep=""), chosen.models)]
      }
      return(assembl.list)
    }
    
    if(em.by == 'all'){
      assembl.list[["TotalConsensus"]] <- chosen.models
      return(assembl.list)
    }
    
    if(em.by == 'PA_dataset+repet'){
      for(dat in .extractModelNamesInfo(chosen.models, info='data.set')){
        for(repet in .extractModelNamesInfo(chosen.models, info='run.eval')){
          mod.tmp <- intersect(x=grep(paste("_",dat,"_",sep=""), chosen.models), 
                               y=grep(paste("_",repet,"_",sep=""), chosen.models))
          if(length(mod.tmp)){
            assembl.list[[paste(dat,"_",repet,'_AllAlgos', sep="")]] <- chosen.models[mod.tmp]
          }
        }
      }
      return(assembl.list)      
    }

    if(em.by == 'PA_dataset+algo'){
      for(dat in .extractModelNamesInfo(chosen.models, info='data.set')){
        for(algo in .extractModelNamesInfo(chosen.models, info='models')){
          mod.tmp <- intersect(x=grep(paste("_",dat,"_",sep=""), chosen.models), 
                               y=grep(paste("_",algo,sep=""), chosen.models))
          if(length(mod.tmp)){
            assembl.list[[paste(dat,"_AllRepet_",algo, sep="")]] <- chosen.models[mod.tmp]
          }
        }
      }
      return(assembl.list)      
    }
  
  }
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  
  em.mod.assemb <- em.models.assembling(chosen.models, em.by)
  
  
  
  ## At this step we made one meta model for each dataset
  ## We have to think about it..
  
#   for(dat in .extractModelNamesInfo(modeling.output@models.computed, info='data.set')){
#     cat("\n   >", dat)

  for(assemb in names(em.mod.assemb) ){
    cat("\n\n  >", assemb, "ensemble modeling")
    models.kept <- em.mod.assemb[[assemb]]
    
    for(eval.m in eval.metric){
      
      if( eval.m != 'none'){
        cat("\n   >", eval.m)
        models.kept.scores <- unlist(lapply(models.kept, function(x){
          mod <- tail(unlist(strsplit(x,"_")), 3)[3]
          run <- tail(unlist(strsplit(x,"_")), 3)[2]
          dat <- tail(unlist(strsplit(x,"_")), 3)[1]
          # select evaluations scores obtained for Evaluation Data if exists or CV if not
          if(modeling.output@has.evaluation.data){
            return(getModelsEvaluations(modeling.output)[eval.m, "Evaluating.data", mod, run, dat])
          } else{
            return(getModelsEvaluations(modeling.output)[eval.m, "Testing.data", mod, run, dat])
          }
          
          }))
        ## set NA to -1
        if(!is.null(models.kept.scores)){
          models.kept.scores[is.na(models.kept.scores)] <- -1
        }
        models.kept <- models.kept[models.kept.scores > eval.metric.quality.threshold[which(eval.metric == eval.m)]]
        models.kept.scores <- models.kept.scores[models.kept.scores > eval.metric.quality.threshold[which(eval.metric == eval.m)]]
      }
      
      if(length(models.kept) ){
        cat("\n   > models kept : ", toString(models.kept))
        if(modeling.output@has.evaluation.data){
          prediction.kept <- as.data.frame(getModelsPredictionEval(modeling.output, as.data.frame = TRUE)[,models.kept])
        } else{
          ## load prediction on each PA dataset
          if(em.by %in% c("PA_dataset",'PA_dataset+algo','PA_dataset+repet')){
            prediction.kept <- as.data.frame(getModelsPrediction(modeling.output, as.data.frame = TRUE)[,models.kept])
          } else{ ## redo prediction on full data.set
            cat("\n   ! Models projection for whole zonation required...")
            prediction.kept <- BIOMOD_Projection(modeling.output = modeling.output,
                                          new.env = getModelsInputData(modeling.output)@data.env.var,
                                          proj.name = 'Tmp',
                                          xy.new.env = getModelsInputData(modeling.output)@coord,
                                          selected.models = models.kept,
                                          compress = 'xz',
                                          clamping.mask = F,
                                          do.stack=T, silent = T)@proj@val
            
            # transform array into data.frame
            prediction.kept <- as.data.frame(prediction.kept)
            names(prediction.kept) <- unlist(lapply(strsplit(names(prediction.kept),".", fixed=TRUE), 
                                                  function(x){
                                                   return(paste(modeling.output@sp.name, x[3], x[2], x[1],sep="_"))
                                                   }))
            # keep only wanted columns
            prediction.kept <- prediction.kept[,models.kept]
            
            prediction.kept <- as.data.frame(prediction.kept)
            unlink(file.path(modeling.output@sp.name,"proj_Tmp"),recursive = TRUE, force = TRUE)
            cat("\n")
          }

        }
      } else {
        cat("\n   ! No models kept due to treshold filtering... Ensemble Modeling was skip!")
        next 
      }

      # 1. Mean of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
      if(prob.mean){
        cat("\n   > Mean of probabilities...")
        em.mean <- round(apply(prediction.kept, 1, mean, na.rm=T))    
      }
      
      # 2. CV of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
      if(prob.cv){
        cat("\n   > Coef of variation of probabilities...")
        em.sd <- apply(prediction.kept, 1, sd, na.rm=T)
        if(!exists('em.mean')) em.mean <- round(apply(prediction.kept, 1, mean, na.rm=T))
        em.cv <- round( em.sd / em.mean,2)
        # putting to 0 points where mean = 0
        em.cv[ em.mean == 0 ] <- 0
      }
      
      # 3. Median of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
      if(prob.median){
        cat("\n   > Median of ptobabilities...")
        em.median <- round(apply(prediction.kept, 1, median, na.rm=T))
      }
      
      # 4. CI of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
      if(prob.ci){
        cat("\n   > Confidence Interval...")
        if(!exists('em.mean')) em.mean <- round(apply(prediction.kept, 1, mean, na.rm=T))
        if(!exists('em.sd')) em.sd <- apply(prediction.kept, 1, sd)
        cat("\n      >", prob.ci.alpha/2*100, "%")
        em.ci.inf <- round(em.mean - qt(1-prob.ci.alpha/2, df = length(models.kept) + 1 ) / sqrt(length(models.kept)) * em.sd)
        cat("\n      >", (1-prob.ci.alpha/2)*100, "%")
        em.ci.sup <- round(em.mean + qt(1-prob.ci.alpha/2, df = length(models.kept) + 1 ) / sqrt(length(models.kept)) * em.sd)
        ### checking pred are on a 0 1000 ladder
        em.ci.inf[em.ci.inf > 1000] <- 1000
        em.ci.inf[em.ci.inf < 0] <- 0
        em.ci.sup[em.ci.sup > 1000] <- 1000
        em.ci.sup[em.ci.sup < 0] <- 0
      }

      # 5. Comitee averaging of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
      if(committee.averaging){      
        cat("\n   >  Comittee averaging...")
        models.kept.tresh <- unlist(lapply(models.kept, function(x){
          mod <- tail(unlist(strsplit(x,"_")), 3)[3]
          run <- tail(unlist(strsplit(x,"_")), 3)[2]
          dat <- tail(unlist(strsplit(x,"_")), 3)[1]
          return(getModelsEvaluations(modeling.output)[eval.m, "Cutoff", mod, run, dat])
          }))
        
        em.ca <- round(apply(as.data.frame(BinaryTransformation(prediction.kept,models.kept.tresh)), 1, mean)*1000)
        ### keep bin thresholds
        EM@em.bin.tresh <- c(EM@em.bin.tresh, 
                             eval(parse( text = paste("list('", modeling.output@sp.name, "_", assemb, "_EMby", eval.m,"' = models.kept.tresh)", sep=""))) )
      }
      
      # 6. weighted mean of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
      if(prob.mean.weight){
        cat("\n   > Prababilities wegthing mean...")
        
        # remove SRE models if ROC
        if(eval.m == 'ROC'){
          sre.id <- grep("_SRE", models.kept)
          if(length(sre.id)>0){
            cat("\n      ! SRE modeling were switched off")
            models.kept <- models.kept[-sre.id]
            models.kept.scores <- models.kept.scores[-sre.id]
            prediction.kept <- prediction.kept[,models.kept]
          }
        }
      
        # weights are "decay" times decreased for each subsequent model in model quality order.                              
        models.kept.scores <- round(models.kept.scores, 10) # sometimes there can be a rounding issue in R, so here I make sure all values are rounded equally.
        
        # dealing with numerical decay
        if(is.numeric(prob.mean.weight.decay)){
          DecayCount <- sum(models.kept.scores>0)
      		WOrder <- order(models.kept.scores, decreasing=T)
      		Dweights <- models.kept.scores
      		for(J in 1:DecayCount) Dweights[WOrder[J]] <- (DecayCount - J + 1) * prob.mean.weight.decay
      		#If 2 or more score are identical -> make a mean weight between the ones concerned
      		for(J in 1:length(models.kept.scores)){
      			if(sum(models.kept.scores[J]==models.kept.scores)>1) Dweights[which(models.kept.scores[J]==models.kept.scores)] <- mean(Dweights[which(models.kept.scores[J]==models.kept.scores)])
      		}      
      		models.kept.scores <- Dweights
      		rm(Dweights,DecayCount,WOrder)          
        }

        ### Standardise model weights
      	models.kept.scores <- models.kept.scores/sum(models.kept.scores)
      	### Compute ensemble forecast
      	em.pmw <- round(as.vector(as.matrix(prediction.kept) %*% models.kept.scores))
        
        ### keep the weighted scores of models
        EM@em.weight <- c(EM@em.weight, eval(parse( text = paste("list('", modeling.output@sp.name, "_", assemb, "_EMby", eval.m,"' = models.kept.scores)", sep=""))))
      }
      
      # 7. Assembling all computed models -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
      em.potential <- c('em.mean','em.cv','em.ci.inf','em.ci.sup','em.median','em.ca','em.pmw')
      em.computed <- em.potential[unlist(lapply(em.potential, exists, envir = environment()))]

      em.pred <- as.data.frame(sapply(em.computed, get, envir = environment()))  
      rm(list=em.computed)
      
      # 8. Evaluating predictions -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
      em.cross.validation <- lapply(em.pred, 
                                    function(x){
                                      if(modeling.output@has.evaluation.data){
                                        obs <-  as.vector(getModelsInputData(modeling.output)@eval.data.species)
                                      } else{
                                        ##### !!!!!! TO DO -> select appropriate part of dataset according to em.by
                                        if(em.by %in% c("PA_dataset",'PA_dataset+algo','PA_dataset+repet')){
                                          if(head(unlist(strsplit(assemb,"_")),1) == 'AllData'){
                                            obs <-  as.vector(getModelsInputData(modeling.output)@data.species)
                                            obs[is.na(obs)] <- 0
                                          } else{
                                            obs <- as.vector(getModelsInputData(modeling.output)@data.species)
                                            obs <- obs[getModelsInputData(modeling.output)@PA[, paste(head(unlist(strsplit(assemb,"_")),1))]]
                                            obs[is.na(obs)] <- 0                                              
                                          }                                
                                        }
                                        
                                        if(em.by %in% c("algo","all") ){
                                          ## we need to take all data even if it is much better to have 
                                          obs <-  as.vector(getModelsInputData(modeling.output)@data.species)
                                          obs[is.na(obs)] <- 0
                                        }                                        

                                      }
#                                         ############################
#                                         cat("\n***")
#                                         return(list(pred = as.vector(x), obs = obs))
#                                         cat("\n***")
#                                     })
#       return(em.cross.validation)
#                                         ############################
                                      em.cr.val <- sapply(unlist(dimnames(getModelsEvaluations(modeling.output))[1]),
                                                      Find.Optim.Stat,
                                                      Fit = as.vector(x),
                                                      Obs = obs,
                                                      Precision = 5)
                                      if(modeling.output@has.evaluation.data){          
                                        rownames(em.cr.val) <- c("Evaluating.data","Cutoff","Sensitivity", "Specificity")
                                      } else{
                                        rownames(em.cr.val) <- c("Testing.data","Cutoff","Sensitivity", "Specificity")
                                      }
                                      em.cr.val <- t(round(em.cr.val,digits=3))
                                      return(em.cr.val)
                                      })
      em.cross.validation <- abind(em.cross.validation, along=3)
      
                                      
      EM@em.computed <- c(EM@em.computed, paste(modeling.output@sp.name, "_", assemb, "_EMby", eval.m,sep=""))
      
      eval(parse(text = paste('EM@em.res$',modeling.output@sp.name, "_", assemb, "_EMby", eval.m, " <- list(em.models.kept = models.kept, em.algo = colnames(em.pred), em.pred = em.pred, em.cross.validation = em.cross.validation )", sep="")))
      
      rm(list=c('em.pred', 'em.cross.validation'))
    } 
  }
  
  .bmCat("Done")  
  return(EM)
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.BIOMOD_EnsembleModeling.check.args <- function(  modeling.output,
                                                   chosen.models,
                                                   eval.metric,
                                                   eval.metric.quality.threshold,
                                                   prob.mean,
                                                   prob.cv,
                                                   prob.ci,
                                                   prob.ci.alpha,
                                                   prob.median,
                                                   committee.averaging,
                                                   prob.mean.weight,
                                                   prob.mean.weight.decay,
                                                   em.by ){
  # 1. modeling.output checking
  if(!(inherits(modeling.output, "BIOMOD.models.out"))){
    stop("Invalid modeling.output argument !\nIt must be a 'BIOMOD.models.out' object")
  }
  
  # 2. chosen.models checking
  if(!length(chosen.models) | (length(chosen.models)==1 & chosen.models[1] == 'all')){ # select all models
    cat("\n   ! all models available will be included in ensemble.modeling")
    chosen.models <- modeling.output@models.computed
  } else{
    chosen.models.check <- chosen.models %in% modeling.output@models.computed
    if(sum(!chosen.models.check) > 0){
      stop(paste("Some selected models not exist : ", toString(chosen.models[!chosen.models.check]),
                 "\nPlease choose models in computed models ( ",
                 toString(modeling.output@models.computed), " )",sep=""))
    }
  }
  
  # 3. eval.metric checking
  if(!is.null(eval.metric)){
    if(!is.character(eval.metric)){
      stop("eval.metric must be a character vector or NULL")
    }
    if('all' %in% eval.metric){
      eval.metric <- dimnames(getModelsEvaluations(modeling.output))[[1]]
    }
    eval.metric.check <- eval.metric %in% dimnames(getModelsEvaluations(modeling.output))[[1]]
    if(sum(!eval.metric.check) > 0){
      stop(paste("Some selected evaluation metrics are not available : ", toString(eval.metric[!eval.metric.check]),
                 "\nPlease choose some in those computed yet ( ",
                 toString(dimnames(getModelsEvaluations(modeling.output))[[1]]), " )",sep=""))
    }
  }
  
  # 4. eval.metric.quality.threshold
  if(!is.null(eval.metric)){
    if(!is.null(eval.metric.quality.threshold)){
      if(!is.numeric(eval.metric.quality.threshold)){
        stop("eval.metric.quality.threshold must be NULL or a numeric vector")
      }
      if(length(eval.metric) != length(eval.metric.quality.threshold)){
        stop("you must give as many eval.metric.quality.threshold than eval.metric (if you give ones)")
      }
      cat("\n   > Evaluation & Weighting methods summary :\n")
      cat(paste(eval.metric, eval.metric.quality.threshold,  sep = " over ", collapse = "\n      "), fill=TRUE, labels = "     ")
    } else{
      cat("\n   ! No eval.metric.quality.threshold -> All models will be kept for Ensemble Modeling")
      eval.metric.quality.threshold <- rep(0, length(eval.metric))
    }
  }
  
  # 5. check selected EM algo
  if( !is.logical(prob.mean) | !is.logical(prob.cv) | !is.logical(prob.ci) | !is.logical(prob.median) |
      !is.logical(committee.averaging) | !is.logical(prob.mean.weight) ){
    stop("prob.mean, prob.cv, prob.ci, prob.median, committee.averaging and prob.mean.weight arguments must be logical")
  }
  if(is.null(eval.metric)){
    if(committee.averaging | prob.mean.weight){
      stop("You must choose eval.metric if you want to compute Comitee Averaging or Probability wegthing mean algorithmes")
    }
  }
  
  # 6. alpha for Confident interval
  if(prob.ci){
    if(!is.numeric(prob.ci.alpha)){
      stop("prob.ci.alpha must be numeric")
    }
    if(prob.ci.alpha <= 0 | prob.ci.alpha>= 0.5){
      stop("prob.ci.alpha must be a numeric between 0 and 0.5")
    }
  }
  
  # 7. decay checking
  if(prob.mean.weight){
    if(is.numeric(prob.mean.weight.decay)){
      if(prob.mean.weight.decay < 0){
        stop("'prob.mean.weight.decay' should be either 'proportional' or a numeric value > 0")
      }
    } else{
      if(prob.mean.weight.decay != 'proportional'){
        stop("'prob.mean.weight.decay' should be either 'proportional' or a numeric value > 0")
      }
#       prob.mean.weight.decay <- 1
    }   
  }

  if(is.null(eval.metric)){
    eval.metric <- 'none'
  }
  
  # 8. by arg checking
  available.em.by <- c('PA_dataset', 'algo', 'all', 'PA_dataset+repet', 'PA_dataset+algo')
  if(!(em.by %in% available.em.by) ){
    stop("Invalid 'em.by' argument given. It must be one of : 'PA_dataset', 'algo', 'all', 'PA_dataset+repet' or 'PA_dataset+algo'")
  }
  
  return( list( modeling.output = modeling.output,
                chosen.models = chosen.models,
                eval.metric = eval.metric,
                eval.metric.quality.threshold = eval.metric.quality.threshold,
                prob.mean = prob.mean,
                prob.cv = prob.cv,
                prob.ci = prob.ci,
                prob.ci.alpha = prob.ci.alpha,
                prob.median = prob.median,
                committee.averaging = committee.averaging,
                prob.mean.weight = prob.mean.weight,
                prob.mean.weight.decay  = prob.mean.weight.decay,
                em.by = em.by))
  
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
                
# .extractModelNamesInfo <- function(model.names, info = 'species'){
#   if(!is.character(model.names)){
#     stop("model.names must be a character vector")
#   }
#   if(!is.character(info) | length(info) != 1 | !(info %in% c('species', 'data.set', 'models', 'run.eval')) ){
#     stop("info must be 'specie', 'data.set', 'models' or 'run.eval'")
#   }
#                 
#   info.tmp <- as.data.frame(strsplit(model.names, "_"))
#   
#   return( switch(info,
#                  species = paste(unique(unlist(info.tmp[-c(nrow(info.tmp), nrow(info.tmp)-1, nrow(info.tmp)-2),])), collapse="_"),
#                  data.set = paste(unique(unlist(info.tmp[(nrow(info.tmp)-2),]))),
#                  run.eval = paste(unique(unlist(info.tmp[(nrow(info.tmp)-1),]))),
#                  models = paste(unique(unlist(info.tmp[(nrow(info.tmp)),]))) ) )
#               
# }