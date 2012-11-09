.Biomod.Models.loop <- function(X,
                                Model,
                                Options,
                                VarImport, 
                                mod.eval.method,
                                SavePred,
                                xy=NULL,
                                rescal.models = TRUE){
  cat("\n\n-=-=-=- Run : ",X$name, '\n')
  res.sp.run <- list()

  for(i in 1:ncol(X$calibLines)){ # loop on RunEval
    cat('\n\n-=-=-=--=-=-=-',paste(X$name,colnames(X$calibLines)[i],sep=""),'\n')
    
    res.sp.run[[colnames(X$calibLines)[i]]] <- lapply(Model, .Biomod.Models, 
                                 Data = X$dataBM,
                                 Options = Options,
                                 calibLines = X$calibLines[,i],
                                 Yweights = X$weights,
                                 nam = paste(X$name,colnames(X$calibLines)[i], sep=""),
                                 VarImport = VarImport,
                                 mod.eval.method = mod.eval.method,
                                 evalData = X$evalDataBM,
                                 SavePred = T,#SavePred,
                                 xy = X$xy,
                                 eval.xy = X$eval.xy,
                                 rescal.models = rescal.models)
    
    names(res.sp.run[[colnames(X$calibLines)[i]]]) <- Model
    
  }
    
  return(res.sp.run)
}


.Biomod.Models <- function (Model, Data, Options, calibLines, Yweights, nam, VarImport = 0, 
                            mod.eval.method = c('ROC','TSS','KAPPA'), evalData = NULL,
                            SavePred = FALSE,
                            xy = NULL, eval.xy = NULL, rescal.models = TRUE){
                              
  ################################################################################################
  # 1. Print model runing and getting model options =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #               
  # check and get modified args if nececary
  timing=FALSE  
  if(timing){
    t0 <- strptime(format(Sys.time(),"%T"), "%T")
  }  
    
  args <- .Biomod.Models.check(Model, Data, Options, calibLines, Yweights, mod.eval.method, evalData)
  
  if(is.null(args)){ # trouble in input data -> Not Run
    return(0)
  } else {
    Data <- args$Data
    Yweights <- args$Yweights
    evalLines <- args$evalLines
    Type <- args$Type
    criteria <- args$criteria
    Prev <- args$Prev
    mod.eval.method <- args$mod.eval.method
    evalData <- args$evalData
  }

  # defining the function outputs
  ListOut <- list(evaluation = NULL,
                  var.import = NULL,
#                   evaluation = ifelse(is.null(mod.eval.method), 
#                                       NA,
#                                       data.frame(matrix(NA, 
#                                                         nrow=length(mod.eval.method), 
#                                                         ncol=4, 
#                                                         dimnames = list(mod.eval.method,
#                                                                         c("Cross.validation","Cutoff",
#                                                                           "Sensitivity", "Specificity"))))),
#                   var.import = ifelse(VarImport<=0,
#                                       NA,
#                                       matrix(0,ncol=(ncol(Data)-1), nrow=1,
#                                              dimnames = list(nam, colnames(Data)[-1]))),
                  pred = NULL,
                  pred.eval = NULL,
                  calib.failure = NULL) 
        
  if (Model == "CTA") {
    
    # converting cost argument
    if(is.null(Options@CTA$cost)){
      cost.tmp <- rep(1,(ncol(Data)-2))
    } else{
      cost.tmp <- Options@CTA$cost
    } 
    if(Options@CTA$parms == 'default'){
      model.sp <- try( rpart(makeFormula(colnames(Data)[1],
                                    head(Data[,-c(1,ncol(Data))]),
                                    'simple', 0),
                        data = Data[calibLines,],
                        weights = Yweights,
                        method = Options@CTA$method,
                        cost = cost.tmp,
                        control = eval(Options@CTA$control)) )    
    } else{
      model.sp <- try( rpart(makeFormula(colnames(Data)[1],
                                    head(Data[,-c(1,ncol(Data))]),
                                    'simple', 0),
                        data = Data[calibLines,],
                        weights = Yweights,
                        method = Options@CTA$method,
                        parms = Options@CTA$parms,
                        cost = cost.tmp,
                        control = eval(Options@CTA$control)) )     
    }
    


    if( !inherits(model.sp,"try-error") ){
      # select best trees --------------- May be done otherway
      tr <- as.data.frame(model.sp$cptable)
      tr$xsum <- tr$xerror + tr$xstd
      tr <- tr[tr$nsplit > 0, ]
      Cp <- tr[tr$xsum == min(tr$xsum), "CP"]
      
      model.sp <- prune(model.sp, cp = Cp[length(Cp)] )
         
      # make prediction
      g.pred <- data.frame(as.integer(as.numeric(predict(model.sp,
                                                         newdata = Data[,-c(1,ncol(Data))],
                                                         type = "prob")[,2]) * 1000))
      
      # rescale or not predictions
      if(rescal.models){
        g.pred <- data.frame(as.integer(.Rescaler5(as.numeric(g.pred[,1]/1000),
                                                   ref = Data[, 1], 
                                                   name = paste(nam,'_',Model,sep=""),
                                                   original = TRUE) * 1000))    
      }
                                                 
      if(!is.null(evalData)){
        # make prediction on evaluation data
        g.pred.eval <- data.frame(as.integer(as.numeric(predict(model.sp, 
                                                         newdata = evalData[,-c(1)],
                                                         type = "prob")[,2]) * 1000))
                                                                
        # rescale or not predictions on evluation data
        if(rescal.models){
          g.pred.eval <- data.frame(as.integer(.Rescaler5(as.numeric(g.pred.eval[,1]/1000),
                                                     ref = evalData[, 1], 
                                                     name = paste(nam,'_',Model,sep=""),
                                                     original = TRUE) * 1000))    
        }                                    
      }

    }
  }
        
        
        
  if (Model == "GAM"){
    
    # NOTE : To be able to take into account GAM options and weights we have to do a eval(parse(...))
    # it's due to GAM implementation ( using of match.call() troubles)
    
    cat("\n\tUser defined control args building..")
    user.control.list <- Options@GAM$control
  
    if(Options@GAM$algo == 'GAM_gam'){
      default.control.list <- gam:::gam.control()
    } else{
      default.control.list <- mgcv:::gam.control()
    }
    
    control.list <- lapply(names(default.control.list), function(x){
      if(x %in% names(user.control.list)){
        return(user.control.list[[x]])
      } else {
        return(default.control.list[[x]])
      }
    })
    names(control.list) <- names(default.control.list)

    ### Old version
    if(Options@GAM$algo == 'GAM_gam'){ ## gam package

      gamStart <- eval(parse(text=paste("gam(",colnames(Data)[1] ,"~1 ," ,
                            " data = Data[calibLines,], family = ", eval(Options@GAM$family),
                            ", weights = Yweights[calibLines])" ,sep="")))
      
      model.sp <- try( step.gam(gamStart, .scope(Data[1:3,-c(1,ncol(Data))], "s", Options@GAM$spline),
                           data = Data[calibLines,],
                           keep = .functionkeep, 
                           direction = "both",
                           trace=control.list$trace,
                           control = eval(control.list)) )

    } else { ## mgcv package
      if(is.null(Options@GLM$myFormula)){
        gam.formula <- makeFormula(colnames(Data)[1],head(Data[,-ncol(Data)]),Options@GAM$type, Options@GAM$interaction.level)
      } else{
        gam.formula <- Options@GLM$myFormula
      }
      
      if (Options@GAM$algo == 'GAM_mgcv'){
        model.sp <- try(mgcv:::gam(gam.formula, 
                            data=Data, 
                            family=Options@GAM$family, 
                            weights = Yweights,
                            control = control.list))
      } else if (Options@GAM$algo == 'BAM_mgcv'){ ## big data.frame gam version
        model.sp <- try(mgcv:::bam(gam.formula, 
                            data=Data, 
                            family=Options@GAM$family,
                            weights = Yweights,
                            control = control.list))
      }
    }

    

    

    if( !inherits(model.sp,"try-error") ){ 
      # make prediction
      g.pred <- data.frame(as.integer(.testnull(model.sp, Prev, Data[,-1]) * 1000))
      
      # rescale or not predictions
      if(rescal.models){
        g.pred <- data.frame(as.integer(.Rescaler5(as.numeric(g.pred[,1]/1000),
                                                   ref = Data[, 1], 
                                                   name = paste(nam,'_',Model,sep=""),
                                                   original = TRUE,
                                                   weights = Yweights) * 1000))    
      }
                                                 
      if(!is.null(evalData)){
        # make prediction on evaluation data
        g.pred.eval <- data.frame(as.integer(.testnull(model.sp, Prev, evalData[,-1]) * 1000))
                                                                
        # rescale or not predictions on evluation data
        if(rescal.models){
          g.pred.eval <- data.frame(as.integer(.Rescaler5(as.numeric(g.pred.eval[,1]/1000),
                                                          ref = evalData[, 1], 
                                                          name = paste(nam,'_',Model,sep=""),
                                                          original = TRUE) * 1000))    
        }                                    
      }
    }
    
  }
        
  if (Model == "GBM") {

      model.sp <- try(gbm(formula = makeFormula(colnames(Data)[1],head(Data)[,-c(1,ncol(Data))], 'simple',0),
                      data = Data[calibLines,], 
                      distribution = Options@GBM$distribution,
                      var.monotone = rep(0, length = ncol(Data)-2), # -2 because of removing of sp and weights
                      weights = Yweights,
                      interaction.depth = Options@GBM$interaction.depth,
                      shrinkage = Options@GBM$shrinkage,
                      bag.fraction = Options@GBM$bag.fraction,
                      train.fraction = Options@GBM$train.fraction,
                      n.trees = Options@GBM$n.trees,
                      verbose = FALSE,
                      cv.folds = Options@GBM$cv.folds))

      if( !inherits(model.sp,"try-error") ){
#     if (exists("model.sp")) {
      best.iter <- try(gbm.perf(model.sp, method = "cv", plot.it = FALSE)) # may be have to be save

      # make prediction
      g.pred <- data.frame(as.integer(predict.gbm(model.sp, Data[,-c(1,ncol(Data))], best.iter,
                                                      type = "response") * 1000))
      
      # rescale or not predictions
      if(rescal.models){
        g.pred <- data.frame(as.integer(.Rescaler5(as.numeric(g.pred[,1]/1000),
                                                   ref = Data[, 1], 
                                                   name = paste(nam,'_',Model,sep=""),
                                                   original = TRUE,
                                                   weights = Yweights) * 1000))    
      }
                                                 
      if(!is.null(evalData)){
        # make prediction on evaluation data
        g.pred.eval <- data.frame(as.integer(predict.gbm(model.sp, evalData[,-c(1)], best.iter,
                                                      type = "response") * 1000))
                                                                
        # rescale or not predictions on evluation data
        if(rescal.models){
          g.pred.eval <- data.frame(as.integer(.Rescaler5(as.numeric(g.pred.eval[,1]/1000),
                                                          ref = evalData[, 1],
                                                          name = paste(nam,'_',Model,sep=""),
                                                          original = TRUE) * 1000))    
        }                                    
      }      
    }
  }
        
        
  if (Model == "GLM"){
    
    ## build the most complete model formula
    if(is.null(Options@GLM$myFormula)){
      glm.formula <- makeFormula(colnames(Data)[1],head(Data),Options@GLM$type, Options@GLM$interaction.level)
    } else{
      glm.formula <- Options@GLM$myFormula
    }
    
    if(Options@GLM$test != 'none'){
      ## make the model selection
      glmStart <- glm(eval(parse(text=paste(colnames(Data)[1],"~1",sep=""))), 
                      data = Data[calibLines,],
                      family = Options@GLM$family,
                      control = eval(Options@GLM$control),
                      weights = Yweights[calibLines],
                      mustart = rep(Options@GLM$mustart, sum(calibLines)),
                      model = TRUE)
      
      ## remove warnings
      warn <- options('warn')                  
      options(warn=-1)
      model.sp <- try( stepAIC(glmStart, 
                          glm.formula,
                          data = Data[calibLines,],
                          direction = "both", trace = FALSE, 
                          k = criteria, 
                          weights = Yweights[calibLines],
                          steps = 10000,
                          mustart = rep(Options@GLM$mustart, sum(calibLines))) ) 
                              
      ## reexec warnings
      options(warn)
                          
    } else {
      ## keep the total model      
      model.sp <- try( glm(glm.formula, 
                      data = cbind(Data[calibLines,],matrix(Yweights[calibLines], ncol=1, dimnames=list(NULL, "Yweights"))) ,
                      family = eval(parse(text=call(Options@GLM$family))),
#                       family = Options@GLM$family,
                      control = eval(Options@GLM$control),
                      weights = Yweights,
#                       mustart = rep(Options@GLM$mustart, sum(calibLines)),
                      model = TRUE) )
    }
                  

                      
    if( !inherits(model.sp,"try-error") ){    
      # print the selected formula
      cat("\n\tselected formula : ")
      print(model.sp$formula, useSource=FALSE)
      
      # make prediction
      g.pred <- data.frame(as.integer(.testnull(model.sp, Prev, Data[,-1]) * 1000))
      
      # rescale or not predictions
      if(rescal.models){
        g.pred <- data.frame(as.integer(.Rescaler5(as.numeric(g.pred[,1]/1000),
                                                   ref = Data[, 1], 
                                                   name = paste(nam,'_',Model,sep=""),
                                                   original = TRUE,
                                                   weights = Yweights) * 1000))    
      }
                                                 
      if(!is.null(evalData)){
        # make prediction on evaluation data
        g.pred.eval <- data.frame(as.integer(.testnull(model.sp, Prev, evalData[,-1]) * 1000))
                                                                
        # rescale or not predictions on evluation data
        if(rescal.models){
          g.pred.eval <- data.frame(as.integer(.Rescaler5(as.numeric(g.pred.eval[,1]/1000),
                                                          ref = evalData[, 1],
                                                          name = paste(nam,'_',Model,sep=""),
                                                          original = TRUE) * 1000))
        }                                    
      }
    }
  
  } 
        
  if (Model == "MARS"){

    model.sp <- try( mars(x = Data[calibLines,2:ncol(Data)],
                     y = Data[calibLines,1],
                     degree = Options@MARS$degree,
                     penalty = Options@MARS$penalty,
                     thresh = Options@MARS$thresh,
                     prune = Options@MARS$prune,
                     w = Yweights[calibLines]) )

    if( !inherits(model.sp,"try-error") ){
      # prediction are automaticly rescaled
      g.pred <- data.frame(as.integer(.Rescaler5(as.numeric(predict(model.sp, Data[,2:ncol(Data)])),
                                                 ref = Data[, 1], 
                                                 name = paste(nam,'_',Model,sep=""),
                                                 original = TRUE,
                                                 weights = NULL#Yweights
                                                 ) * 1000))
                                                 
      if(!is.null(evalData)){
        g.pred.eval <- data.frame(as.integer(.Rescaler5(as.numeric(predict(model.sp, evalData[,2:ncol(evalData)])),
                                                 ref = evalData[, 1], 
                                                 name = paste(nam,'_',Model,sep=""),
                                                 original = TRUE) * 1000))
      }
    }
  }
        
  if (Model == "FDA") {
    model.sp <- try(fda(formula = makeFormula(colnames(Data)[1],head(Data)[,-c(1,ncol(Data))], 'simple',0),
                    data = Data[calibLines,],
                    method = eval(parse(text=call(Options@FDA$method))),
                    weights = Yweights))
                    
    if( !inherits(model.sp,"try-error") ){
      # prediction are automaticly rescaled
      g.pred <- try(data.frame(as.integer(.Rescaler5(as.numeric(predict(model.sp, Data[,-c(1,ncol(Data))],
                                                                        type = "posterior")[, 2]),
                                                     ref = as.numeric(Data[, 1]),
                                                     name = paste(nam,'_',Model,sep=""),
                                                     original = TRUE,
                                                     weights = Yweights) * 1000)))
      
      if( inherits(g.pred,"try-error") ){
        # remove g.pred if not run
        rm('g.pred')
      }
      
      if(!is.null(evalData)){
        g.pred.eval <- data.frame(as.integer(.Rescaler5(as.numeric(predict(model.sp, evalData[,-c(1)],
                                                 type = "posterior")[, 2]),
                                                 ref = as.numeric(evalData[, 1]), 
                                                 name = paste(nam,'_',Model,sep=""),
                                                 original = TRUE) * 1000))
      }
    }
  }
            
  if (Model == "ANN") {
      CV_nnet = .CV.nnet(Input = Data[,2:(ncol(Data)-1)], 
                         Target = Data[calibLines,1], 
                         nbCV = Options@ANN$NbCV, 
                         W = Yweights[calibLines])
      
      model.sp <- try(nnet(formula = makeFormula(colnames(Data)[1],head(Data[,-c(1,ncol(Data))]), 'simple',0),
                       data = Data[calibLines,],
                       size = CV_nnet[1,1],
                       rang = Options@ANN$rang,
                       decay = CV_nnet[1, 2],
                       weights=Yweights,
                       maxit = Options@ANN$maxit,
                       trace = FALSE))

      if( !inherits(model.sp,"try-error") ){
        # prediction are automaticly rescaled
        g.pred <- data.frame(as.integer(.Rescaler5(as.numeric(predict(model.sp, Data[,-c(1,ncol(Data))], type = "raw")),
                                                   ref = Data[,1],
                                                   name = paste(nam,'_',Model,sep=""),
                                                   original = TRUE,
                                                   weights = Yweights) * 1000))
                                                     
        if(!is.null(evalData)){
         g.pred.eval <- data.frame(as.integer(.Rescaler5(as.numeric(predict(model.sp, evalData[,-c(1)], type = "raw")),
                                                   ref = evalData[,1], 
                                                   name = paste(nam,'_',Model,sep=""),
                                                   original = TRUE) * 1000))
        }
      }
  }

        
  if (Model == "RF") {
    
    if(Options@RF$do.classif){
      # defining occurences as factor for doing classification and not regression in RF
      dTmp <- Data[,1]
      Data[,1] <- as.factor(Data[,1])  
    }
    
    if(Options@RF$mtry == 'default'){
      model.sp <- try(randomForest(formula = makeFormula(colnames(Data)[1],head(Data), 'simple',0),
                               data = Data[calibLines,],
                               ntree = Options@RF$ntree,
                               #mtry = ifelse(Options@RF$ntree == 'default', round((ncol(Data)-1)/2), Options@RF$ntree ),
                               importance = FALSE,
                               norm.votes = TRUE,
                               strata = factor(c(0,1))) )      
    } else {
      model.sp <- try(randomForest(formula = makeFormula(colnames(Data)[1],head(Data), 'simple',0),
                               data = Data[calibLines,],
                               ntree = Options@RF$ntree,
                               mtry = Options@RF$mtry,
                               importance = FALSE,
                               norm.votes = TRUE,
                               strata = factor(c(0,1))) )
    }


    if(Options@RF$do.classif){                             
      # canceling occurences class modifications
      Data[,1] <- as.numeric(dTmp)
      rm(dTmp)
    }
    
    if( !inherits(model.sp,"try-error") ){
      # make prediction
      g.pred <- data.frame(as.integer(as.numeric(predict(model.sp,Data[,-1], type='prob')[,'1']) *1000))
      
      # rescale or not predictions
      if(rescal.models){
        g.pred <- data.frame(as.integer(.Rescaler5(as.numeric(g.pred[,1]/1000),
                                                   ref = Data[, 1], 
                                                   name = paste(nam,'_',Model,sep=""),
                                                   original = TRUE,
                                                   weights = Yweights) * 1000))    
      }
                                                 
      if(!is.null(evalData)){
        # make prediction on evaluation data
        g.pred.eval <- data.frame(as.integer(as.numeric(predict(model.sp,evalData[,-1], type='prob')[,'1']) *1000))
                                                                
        # rescale or not predictions on evluation data
        if(rescal.models){
          g.pred.eval <- data.frame(as.integer(.Rescaler5(as.numeric(g.pred.eval[,1]/1000),
                                                     ref = evalData[, 1], 
                                                     name = paste(nam,'_',Model,sep=""),
                                                     original = TRUE) * 1000))
        }                                    
      }
    }
    
  }

  if (Model == "SRE"){
    g.pred <- data.frame(as.integer(as.numeric(sre(Data[calibLines,1],
                                                   Data[calibLines,2:ncol(Data)],
                                                   Data[,2:ncol(Data)], 
                                                   Options@SRE$quant)) * 1000))
    if(!is.null(evalData)){
      g.pred.eval <- data.frame(as.integer(as.numeric(sre(Data[calibLines,1],
                                                   Data[calibLines,2:ncol(Data)],
                                                   evalData[,2:ncol(evalData)], 
                                                   Options@SRE$quant)) * 1000))
    }

    # save data and quant Used for this prediction
    dir.create(paste(getwd(), .Platform$file.sep, colnames(Data)[1], .Platform$file.sep, "models", .Platform$file.sep, nam, "_SRE", sep=""), showWarnings=F)
    eval(parse(text=paste("Data_", nam, "_SRE = list()
                          Data_", nam, "_SRE$Response <- ",
                          "Data[calibLines,1]", "
                          Data_", nam, "_SRE$Explanatory <- ",
                          "Data[calibLines,2:ncol(Data)]","
                          Data_", nam, "_SRE$Quant <- ",
                          Options@SRE$quant,"
                          save(Data_", nam, "_SRE, file='",
                          getwd(), .Platform$file.sep, colnames(Data)[1], .Platform$file.sep, "models", .Platform$file.sep, nam, "_SRE", .Platform$file.sep, "Data_", nam, "_SRE')",
                          sep="")))

  }
                       
  if (Model == "MAXENT"){
    .Prepare.Maxent.WorkDir(Data, xy, calibLines, nam, VarImport, evalData, eval.xy, species.name=colnames(Data)[1])
    cat("\n Runing Maxent...")

    # run MaxEnt:
#     if(VarImport > 0 ){ # projection on suffle table to be able to compute VarImport latter
#       if(.Platform$OS.type == "windows"){
#         currentDir <- gsub(pattern=" ", replacement="~1", x=getwd())
#       } else{
#         currentDir <- getwd()
#       }
    
      system(command=paste("java -mx512m -jar ", file.path(Options@MAXENT$path_to_maxent.jar, "maxent.jar"), " environmentallayers=\"",
                           file.path(getwd(), colnames(Data)[1], "MaxentTmpData", "Back_swd.csv"),"\" samplesfile=\"",
                           file.path(getwd(), colnames(Data)[1], "MaxentTmpData", "Sp_swd.csv"),"\" projectionlayers=\"",
                           gsub(", ",",",toString(list.files(paste(getwd(), .Platform$file.sep, colnames(Data)[1], .Platform$file.sep, "MaxentTmpData", .Platform$file.sep, "Pred",sep=""),
                                               full.names= T))), "\" outputdirectory=\"",
                           file.path(getwd(), colnames(Data)[1], "models", paste(nam, "_MAXENT_outputs", sep="")),"\"",
                           " outputformat=logistic ",
#                            "jackknife maximumiterations=",Options@MAXENT$maximumiterations, 
                           " redoifexists",
                           " visible=", Options@MAXENT$visible,
                           " linear=", Options@MAXENT$linear,
                           " quadratic=", Options@MAXENT$quadratic,
                           " product=", Options@MAXENT$product,
                           " threshold=", Options@MAXENT$threshold,
                           " hinge=", Options@MAXENT$hinge,
                           " lq2lqptthreshold=", Options@MAXENT$lq2lqptthreshold,
                           " l2lqthreshold=", Options@MAXENT$l2lqthreshold,
                           " hingethreshold=", Options@MAXENT$hingethreshold,
                           " beta_threshold=", Options@MAXENT$beta_threshold,
                           " beta_categorical=", Options@MAXENT$beta_categorical,
                           " beta_lqp=", Options@MAXENT$beta_lqp,
                           " beta_hinge=", Options@MAXENT$beta_hinge,
                           " defaultprevalence=", Options@MAXENT$defaultprevalence,
                           " autorun nowarnings notooltips", sep=""), wait = TRUE)
    
#     } else{ # classic Run
#       system(command=paste("java -mx4000m -jar maxent.jar environmentallayers=",
#                            getwd(),"/MaxentTmpData/Back_swd.csv samplesfile=",
#                            getwd(),"/MaxentTmpData/Sp_swd.csv outputdirectory=",
#                            getwd(),"/",colnames(Data)[1],"/models/", nam, "_MAXENT outputformat='logistic' ",
# #                            "jackknife maximumiterations=",Options@MAXENT$maximumiterations, 
#                            " redoifexists ",
#                            "autorun nowarnings notooltips", sep=""))      
#     }
    
#     g.pred <- read.csv(paste(getwd(),"/",colnames(Data)[1], "/models/", nam, "_MAXENT/", 
#                                nam,".csv", sep=""))[,3]
    
    # create a maxent_model object
    model.sp <- maxent_model(sp_name =  colnames(Data)[1], 
                             model_name=paste(nam,'_',Model,sep=""), 
                             model_output_dir=file.path(getwd(), colnames(Data)[1], "models", paste(nam, "_MAXENT_outputs", sep="")),
                             model_options=Options@MAXENT)
    
    cat("\n Getting predictions...")
    g.pred <- read.csv(file.path(getwd(), colnames(Data)[1], "models", paste(nam,'_MAXENT_outputs',sep=""), paste(nam,"_Pred_swd.csv", sep="") ) )[,3]
                   
    g.pred <- data.frame(as.integer(.Rescaler5(g.pred,
                                               ref = Data[,1], 
                                               name = paste(nam,'_',Model,sep=""),
                                               original = TRUE,
                                               weights = Yweights) * 1000))

    if(!is.null(evalData)){
      g.pred.eval <- read.csv(paste(getwd(), .Platform$file.sep, colnames(Data)[1], .Platform$file.sep, "models", .Platform$file.sep,  nam, "_MAXENT_outputs", .Platform$file.sep, nam,"_Pred_eval_swd.csv", sep=""))[,3]
                   
      g.pred.eval <- data.frame(as.integer(.Rescaler5(g.pred.eval,
                                               ref = evalData[,1], 
                                               name = paste(nam,'eval_',Model,sep=""),
                                               original = TRUE) * 1000))
                                                      
      file.remove(list.files(path=paste(getwd(), .Platform$file.sep, colnames(Data)[1], .Platform$file.sep,
                                                  "models", .Platform$file.sep, nam, "_MAXENT_outputs", .Platform$file.sep, sep=""),
                                       pattern=paste(nam,"_Pred_eval_swd", sep=""),
                                       full.names=TRUE))
    }
         
  }

  if(timing){
    t1 <- strptime(format(Sys.time(),"%T"), "%T")
  }
    
  if (!exists("g.pred")) { # keep the name of uncompleted modelisations
    cat("\n   ! Note : ",paste(nam,'_',Model,sep=''), "failed!\n")
    ListOut$calib.failure = paste(nam,'_',Model,sep='')
    return(ListOut)
  }

  if (exists("g.pred")){

    
    ## Evaluation result stuff
    if(length(mod.eval.method) > 0){
      cat("\n\tEvaluating Model stuff...")
    }
    
    ## Check no NA in g.pred to avoid evaluation failures
    if(sum(is.na(g.pred)) > 0 ){
      g.pred.without.na <- data.frame(g.pred[-which(is.na(g.pred)),1])
      cat('\n\tNote : some NA occurs in predictions')
    } else {
      g.pred.without.na <- data.frame(g.pred[,1])
    }
    
    cross.validation <- sapply(mod.eval.method,
                              Find.Optim.Stat,
                              Fit = g.pred.without.na[evalLines,],
                              Obs = Data[evalLines,1],
                              Precision = 5)
    rownames(cross.validation) <- c("Testing.data","Cutoff","Sensitivity", "Specificity")
    
    if(exists('g.pred.eval')){
      ## Check no NA in g.pred to avoid evaluation failures
      if(sum(is.na(g.pred)) > 0 ){
        g.pred.eval.without.na <- data.frame(g.pred.eval[-which(is.na(g.pred.eval)),1])
        cat('\n\tNote : some NA occurs in predictions')
      } else {
        g.pred.eval.without.na <- data.frame(g.pred.eval[,1])
      }
      
      true.evaluation <- sapply(mod.eval.method,
                                function(x){
                                  return( Find.Optim.Stat(Stat = x,
                                                          Fit = g.pred.eval.without.na[,],
                                                          Obs = evalData[,1],
                                                          Fixed.thresh = cross.validation["Cutoff",x]) )
                                })
      
#       rownames(true.evaluation) <- c("Evaluating.data","Cutoff","Sensitivity", "Specificity")
      
      cross.validation <- rbind(cross.validation["Testing.data",], true.evaluation)
      
      rownames(cross.validation) <- c("Testing.data","Evaluating.data","Cutoff","Sensitivity", "Specificity")
    }
    
    ListOut$evaluation <- t(round(cross.validation,digits=3))
    
    rm(list=c('cross.validation', 'g.pred.without.na') )
                               
  } # end evaluation stuff 
  
  if(timing){
    t2 <- strptime(format(Sys.time(),"%T"), "%T")
  }    

  if(exists("g.pred")){ #saving Predictions
    if(SavePred) ListOut$pred <- g.pred[, ]

    # keep the model name
    ListOut$ModelName <- paste( nam, "_", Model, sep = "")
  }

  if(exists("g.pred.eval")){ # saving Predictions on evaluation data
    ListOut$pred.eval <- g.pred.eval[, ]
  }
                                  
  if (!(Model %in% c("SRE"))){ # saving Models
   eval( parse( text = paste( nam, "_", Model, " = model.sp ", sep = "")))       
    
   eval( parse( text = paste("save(",nam, "_", Model, ",file='", getwd(), .Platform$file.sep,
                             unlist(strsplit(nam,'_'))[1] , .Platform$file.sep, "models", .Platform$file.sep,
                             nam, "_", Model, "', compress='",
                             ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz'),
                             "')", sep = "")))
   if(Model == 'GBM'){
     if(exists('best.iter'))
       save(best.iter, file=paste( getwd(), .Platform$file.sep, unlist(strsplit(nam,'_'))[1] , .Platform$file.sep, "models",  .Platform$file.sep, nam, "_", Model,"_best.iter", sep=""), compress=ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz'))
   }
   
   
      
  } #  end saving Models


     if (VarImport > 0){ # do Varimp stuff
       
       # Create data.frame vhere corelation between predictions with and without permutation will be 
       # stored
       
        cat("\n\tEvaluating Predictor Contributions...", "\n")
        TempVarImp <- as.data.frame(matrix(data = 0, nrow = 1,ncol = (ncol(Data) - 1 )) )
        colnames(TempVarImp) <- colnames(Data)[-1]

        if (Model %in% c("CTA","ANN","GBM","FDA")){TempVarImp <- TempVarImp[,-ncol(TempVarImp)]} # remove weights column
     
        for (J in 1:(ncol(Data) - 1 )) {
            for (K in 1:VarImport) {
              TempDS <- Data[, -1]
              TempDS[, J] <- sample(TempDS[, J])
              
                if (Model == "ANN") {
                  if(J < ncol(TempDS)){
                    set.seed(555)
                    TempVarImp[1, J] <- TempVarImp[1, J] + cor(g.pred[, 
                      ], as.integer(.Rescaler5(as.numeric(predict(model.sp, 
                      TempDS, type = "raw")), ref = Data[, 1], 
                                              name = paste(nam,'_',Model,sep="")) * 1000))
                  }
                }
  
                if (Model == "CTA"){
                  if(J < ncol(TempDS)){
                    TempVarImp[1, J] <- TempVarImp[1, J] + cor(g.pred[, 
                      ], as.integer(as.numeric(predict(model.sp, 
                      TempDS, type='prob')[,2] * 1000)))
                  } 
                }
              
                if (Model %in% c("GLM","GAM") )
                  TempVarImp[1, J] <- TempVarImp[1, J] + cor(g.pred[, 
                    ], as.integer(as.numeric(.testnull(model.sp, 
                    Prev, TempDS)) * 1000))

                if (Model == "GBM"){
                  if(J < ncol(TempDS)){
                    TempVarImp[1, J] <- TempVarImp[1, J] + cor(g.pred[, 
                      ], as.integer(as.numeric(predict.gbm(model.sp, 
                      TempDS, best.iter, type = "response")) * 1000))
                  }
                }
              
                if (Model == "MARS")
                  TempVarImp[1, J] <- TempVarImp[1, J] + cor(g.pred[, 
                    ], as.integer(.Rescaler5(predict(model.sp, 
                    TempDS), ref = Data[,1], name = paste(nam,'_',Model,sep="")) * 1000))
              
                if (Model == "FDA"){ 
                  if(J < ncol(TempDS)){
                    TempVarImp[1, J] <- TempVarImp[1, J] + cor(g.pred[, 
                      ], as.integer(.Rescaler5(predict(model.sp, TempDS, type = "posterior")[, 2],
                                               ref = Data[,1],
                                               name = paste(nam,'_',Model,sep="")) * 1000))
                  }                                             
                }
              
                if (Model == "RF") 
                  TempVarImp[1, J] <- TempVarImp[1, J] + cor(g.pred[, 
                    ], as.integer(as.numeric(predict(model.sp, 
                    TempDS)) * 1000))
             
                if (Model == "SRE") 
                  TempVarImp[1, J] <- TempVarImp[1, J] + cor(g.pred[, 
                    ], as.integer(as.numeric(sre(Data[,1], TempDS, 
                    Data[,2:ncol(Data)], Options@SRE$quant)) * 1000))
              
              if(Model == "MAXENT"){
                # we have created all the permutation at model building step
                var.imp.tmp <- read.csv(paste(getwd(),.Platform$file.sep,colnames(Data)[1], .Platform$file.sep,
                                              "models", .Platform$file.sep, nam, "_MAXENT_outputs", .Platform$file.sep,
                                              nam, "_", colnames(Data)[J+1], "_", K, "_swd.csv", sep=""))

                TempVarImp[1, J] <- TempVarImp[1, J] + cor(g.pred[,],
                                                           as.integer(.Rescaler5(var.imp.tmp[,3],
                                                                      ref = Data[,1],
                                                                      name = paste(nam,'_',Model,sep="")) * 1000))
                # removing useless files from workspace and hardDisk
                rm(var.imp.tmp)
                file.remove(list.files(path=paste(getwd(), .Platform$file.sep, colnames(Data)[1], .Platform$file.sep,
                                                  "models", .Platform$file.sep, nam, "_MAXENT_outputs", .Platform$file.sep, sep=""),
                                       pattern=paste(nam,"_",colnames(Data)[J+1],"_",K,"_swd", sep=""),
                                       full.names=TRUE))
              }
                  
            }
        }
                                                               
        ListOut$var.import <- sapply(round(as.numeric(1 - (TempVarImp/VarImport)), digits = 3),min,1)

     } # end do Varimp stuff
  
    if(timing){
      t3 <- strptime(format(Sys.time(),"%T"), "%T")
      try(cat(Model,"\t",difftime(t1,t0, units='secs')[[1]],"\t",difftime(t2,t1, units='secs')[[1]],"\t",difftime(t3,t2, units='secs')[[1]],'\n', file="Biomod_Modelling_Timing.txt", append=T ))
    }
     
  return(ListOut)
}


.Biomod.Models.check <- function(Model, Data, Options, calibLines, Yweights, mod.eval.method, evalData, criteria=NULL, Prev=NULL){
  # replace Pseudo absences selected (NA) into true absences (0).. for model computing purpose

  if(sum(is.na(Data[,1])))
    Data[which(is.na(Data[,1])),1] <- 0
  
  # Calib Lines checking
  # & Test if there is absences AND presences in data given
  if (sum(!calibLines)>0){ # data are splited into 2 set one for evaluation and an other for evaluation stuff
    evalLines <- !calibLines
    if(sum(Data[calibLines,1] == 0 ) == 0 || sum(Data[calibLines,1] == 0 ) == sum(calibLines) ||
      sum(Data[evalLines,1] == 0) == 0 || sum(Data[evalLines,1] == 0) == sum(evalLines)){
     warning(paste(colnames(Data)[1], " ", Model," was switch off because of no both
                   presences and absences data given",sep=""), immediate.=T)
     return(NULL)
    }
  } else { # all values are taken for cali and valid stuff -----> Not so good but usefull for little data set
    evalLines <- calibLines
    if(sum(Data[,1] == 0 ) == 0 || sum(Data[,1] == 0 ) == nrow(Data)){
      warning(paste(colnames(Data)[1], " ", Model," was switch off because of no both
                   presences and absences data given (full model)",sep=""), immediate.=T)
      return(NULL)
    }
  }
  
  # weights checking
  if(is.null(Yweights)){
    Yweights <- rep(1,nrow(Data))
  }
  
  if(Model %in% c('GBM','CTA','ANN','FDA','GAM')){ # this models required data and weights to be in a same datdaset
    Data <- cbind(Data,Yweights)
  }
  
  # models options checking and printing
  if (Model == "GLM"){
    cat('\nModel=GLM')
    if(!is.null(Options@GLM$myFormula)){
      cat('\n\tformula = ', paste(Options@GLM$myFormula[2],Options@GLM$myFormula[1],Options@GLM$myFormula[3]))
    } else{
      cat('',Options@GLM$type,'with', ifelse(Options@GLM$interaction.level == 0, 'no interaction', paste('order',Options@GLM$interaction.level,'interaction level')))                  
    }

    if(Options@GLM$test == "AIC"){
      criteria <- 2
      cat("\n\tStepwise procedure using AIC criteria")      
    } else if(Options@GLM$test == "BIC"){
      criteria <- log(ncol(Data))
      cat("\n\tStepwise procedure using BIC criteria")     
    } else {#if(Options@GLM$test == "none"){
      cat("\n\tNo stepwise procedure")
    }
      
  }
                 
    if (Model == "GBM") {
        cat("\nModel=Generalised Boosting Regression \n")
        cat("\t", Options@GBM$n.trees, "maximum different trees and ", Options@GBM$cv.folds,
            " Fold Cross-Validation")
        set.seed(456) # to be able to refind our trees MAY BE BAD
    }
 
    if (Model == "GAM") {
        cat("\nModel=GAM")
        cat("\n\t",Options@GAM$algo,"algorithm chosen")
#         cat("\t", Options@GAM$spline, " Degrees of smoothing")
#         if(is.null(Yweights)) Yweights <- rep(1,nrow(Data))
#         Data <- cbind(Data,Yweights)
    }
 
    if (Model == "CTA") {
        cat("\nModel=Classification tree \n")
        cat("\t", Options@CTA$control$xval, "Fold Cross-Validation")
        set.seed(123) # to be able to refind our trees MAY BE BAD
    }
 
    if (Model == "ANN") {
        cat("\nModel=Artificial Neural Network \n")
        cat("\t", Options@ANN$NbCV, "Fold Cross Validation + 3 Repetitions")
#         cat("\tCalibration and evaluation phase: Nb of cross-validations: ", 
#             ncol(Ids), "\n")
        set.seed(555) # to be able to refind our trees MAY BE BAD
    }

    if (Model == "SRE") 
        cat("\nModel=Surface Range Envelop")

    if (Model == "FDA"){ 
        cat("\nModel=Flexible Discriminant Analysis")
    }
  
    if (Model == "MARS"){ 
        cat("\nModel=Multiple Adaptive Regression Splines")
    }
  
    if (Model == "RF"){ 
        cat("\nModel=Breiman and Cutler's random forests for classification and regression")
        set.seed(71)
    }
  
    if(Model == 'MAXENT'){
      cat('\nModel=MAXENT')
    }
  
#     if (Model == "GLM" | Model == "GAM") 
#         Prev <- sum(DataBIOMOD[, i + Biomod.material$NbVar])/nrow(DataBIOMOD)
    ## not exactly same as before
  if (Model == "GLM" | Model == "GAM"){ 
    Prev <- sum(Data[,1])/length(Data[,1])
  }
  
  # Evaluation Check
  available.eval.meth <- c('ROC','KAPPA','TSS','ACCURACY','BIAS','POD','FAR','POFD','SR','CSI',
                           'ETS','HK','HSS','OR','ORSS')
  
#   if( Model %in% c('SRE') ) available.eval.meth <- available.eval.meth[which(available.eval.meth!='ROC')]
  if(sum(!(mod.eval.method %in% available.eval.meth)) > 0 ){
    warnings(paste(toString(mod.eval.method[!which(mod.eval.method %in% available.eval.meth)]),
                   ' were switched off !', sep='' ), imediate = TRUE)
  }
  mod.eval.method <- mod.eval.method[which(mod.eval.method %in% available.eval.meth)]
   
  return(list(Data=Data,
              Yweights=Yweights,
              evalLines=evalLines,
              criteria=criteria,
              Prev=Prev,
              mod.eval.method=mod.eval.method,
              evalData=evalData))

}