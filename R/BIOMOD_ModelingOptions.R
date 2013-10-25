####################################################################################################
# BIOMOD_ModelingOptions
# Damien.G
# feb 2012
####################################################################################################

# AIM :
#   Define Models options if defaut configuration is not appropreated

# INPUT :


# OUTPUT : 



# NOTE : 


####################################################################################################
'BIOMOD_ModelingOptions' <- function(
                        GLM = NULL, 
                        GBM = NULL,
                        GAM = NULL,
                        CTA = NULL,
                        ANN = NULL,
                        SRE = NULL,
                        FDA = NULL,
                        MARS = NULL,
                        RF = NULL,
                        MAXENT = NULL
                        ){
  # 1. create a defaut BIOMOD.Model.Options object
  opt <- new('BIOMOD.Model.Options')
  
  # 2. modify it if necessary
  if(!is.null(GLM)){
    if(!is.null(GLM$type)) { opt@GLM$type <- GLM$type }
    if(!is.null(GLM$interaction.level)) { opt@GLM$interaction.level <- GLM$interaction.level }
    if(!is.null(GLM$myFormula)) { opt@GLM$myFormula <- GLM$myFormula }
    if(!is.null(GLM$test)) { opt@GLM$test <- GLM$test }
    if(!is.null(GLM$family)) {
      fam.test <- TRUE
      if( class(GLM$family) == 'family'){
        opt@GLM$family <- GLM$family
      } else{
        if( is.character(GLM$family)){
          if(! unlist(strsplit(GLM$family,"[/(]"))[1] %in% c('binomial', 'gaussian', 'Gamma', 'inverse.gaussian', 'poisson', 'quasi', 'quasibinomial', 'quasipoisson')){ fam.test <- FALSE}
          
          if(grepl(')', GLM$family)){ # check string formalisation to add () if necessary
            opt@GLM$family <- eval(parse(text=GLM$family))
          } else{
            opt@GLM$family <- eval(parse(text=paste(GLM$family,"()", sep="")))
          }
        } else{ fam.test <- FALSE }
      }
      if(!fam.test){
        cat("\n!!! invalid GLM$family given -> binomial(link = 'logit') was automatically set up !!!")
        opt@GLM$family <- binomial(link = 'logit')
      }
    }
    if(!is.null(GLM$mustart)) { opt@GLM$mustart <- GLM$mustart }
    if(!is.null(GLM$control)) { opt@GLM$control <- GLM$control }
  }
  
  if(!is.null(GBM)){
#     if(!is.null(GBM$type )) { opt@GBM$type <- GBM$type }
#     if(!is.null(GBM$interaction.level )) { opt@GBM$interaction.level <- GBM$interaction.level }
    if(!is.null(GBM$distribution )) { opt@GBM$distribution <- GBM$distribution }
    if(!is.null(GBM$n.trees )) { opt@GBM$n.trees <- GBM$n.trees }
    if(!is.null(GBM$interaction.depth )) { opt@GBM$interaction.depth <- GBM$interaction.depth }
    if(!is.null(GBM$n.minobsinnode )) { opt@GBM$n.minobsinnode <- GBM$n.minobsinnode }
    if(!is.null(GBM$shrinkage )) { opt@GBM$shrinkage <- GBM$shrinkage }
    if(!is.null(GBM$bag.fraction )) { opt@GBM$bag.fraction <- GBM$bag.fraction }
    if(!is.null(GBM$train.fraction )) { opt@GBM$train.fraction <- GBM$train.fraction }
    if(!is.null(GBM$cv.folds )) { opt@GBM$cv.folds <- GBM$cv.folds }
    if(!is.null(GBM$keep.data )) { opt@GBM$keep.data <- GBM$keep.data }
    if(!is.null(GBM$verbose )) { opt@GBM$verbose <- GBM$verbose }
#     if(!is.null(GBM$class.stratify.cv )) { opt@GBM$class.stratify.cv <- GBM$cv.folds }
    if(!is.null(GBM$perf.method )) { opt@GBM$perf.method <- GBM$perf.method }
  }


  
  if(!is.null(GAM)){
    if(!is.null(GAM$algo )) { opt@GAM$algo <- GAM$algo }
    if(!is.null(GAM$type )) { opt@GAM$type <- GAM$type }
    if(!is.null(GAM$k )) { opt@GAM$k <- GAM$k } else{
      if(opt@GAM$algo == 'GAM_gam'){
        opt@GAM$k <- 4
      } else{
        opt@GAM$k <- -1
      }
    }
    if(!is.null(GAM$interaction.level )) { opt@GAM$interaction.level <- GAM$interaction.level }
    if(!is.null(GAM$myFormula )) { opt@GAM$myFormula <- GAM$myFormula }
    if(!is.null(GAM$family)) {
      fam.test <- TRUE
      if( class(GAM$family) == 'family'){
        opt@GAM$family <- GAM$family
      } else{
        if( is.character(GAM$family)){
          if(! unlist(strsplit(GAM$family,"[/(]"))[1] %in% c('binomial', 'gaussian', 'Gamma', 'inverse.gaussian', 'poisson', 'quasi', 'quasibinomial', 'quasipoisson')){ fam.test <- FALSE}
          
          if(grepl(')', GAM$family)){ # check string formalisation to add () if necessary
            opt@GAM$family <- eval(parse(text=GAM$family))
          } else{
            opt@GAM$family <- eval(parse(text=paste(GAM$family,"()", sep="")))
          }
        } else{ fam.test <- FALSE }
      }
      if(!fam.test){
        cat("\n!!! invalid GAM$family given -> binomial(link = 'logit') was automatically set up !!!")
        opt@GAM$family <- binomial(link = 'logit')
      }
    }
    
    if(is.null(GAM$control )) {
      if(opt@GAM$algo == 'GAM_gam'){
        opt@GAM$control <- gam::gam.control()
      } else{ opt@GAM$control <- mgcv::gam.control() }
    } else{
      user.control.list <- GAM$control
      if(opt@GAM$algo == 'GAM_gam'){
        default.control.list <- gam::gam.control()
      } else{
        default.control.list <- mgcv::gam.control()
      }
      
      control.list <- lapply(names(default.control.list), function(x){
        if(x %in% names(user.control.list)){
          return(user.control.list[[x]])
        } else {
          return(default.control.list[[x]])
        }
      })
      
      names(control.list) <- names(default.control.list)
      opt@GAM$control <- control.list
    }
      
    if(!is.null(GAM$method )) { opt@GAM$method <- GAM$method }
    if(!is.null(GAM$optimizer )) { opt@GAM$optimizer <- GAM$optimizer }
    if(!is.null(GAM$select )) { opt@GAM$select <- GAM$select }
    if(!is.null(GAM$knots )) { opt@GAM$knots <- GAM$knots }
    if(!is.null(GAM$paraPen )) { opt@GAM$paraPen <- GAM$paraPen } 
  } else{
    if(opt@GAM$algo == 'GAM_gam'){
      opt@GAM$control <- gam::gam.control()
      opt@GAM$k <- 4
    } else{ 
      opt@GAM$control <- mgcv::gam.control()
      opt@GAM$k <- -1
    }
  }

  
  if(!is.null(CTA)){
#     if(!is.null(CTA$type )) { opt@CTA$type <- CTA$type }
#     if(!is.null(CTA$interaction.level )) { opt@CTA$interaction.level <- CTA$interaction.level }
    if(!is.null(CTA$method )) { opt@CTA$method <- CTA$method }
    if(!is.null(CTA$parms )) { opt@CTA$parms <- CTA$parms }
    if(!is.null(CTA$control )) { opt@CTA$control <- CTA$control }
    if(!is.null(CTA$cost )) { opt@CTA$cost <- CTA$cost }
  }
   
  if(!is.null(ANN)){
#     if(!is.null(ANN$type )) { opt@ANN$type <- ANN$type }
#     if(!is.null(ANN$interaction.level )) { opt@ANN$interaction.level <- ANN$interaction.level }
    if(!is.null(ANN$NbCV )) { opt@ANN$NbCV <- ANN$NbCV }
    if(!is.null(ANN$rang )) { opt@ANN$rang <- ANN$rang }
    if(!is.null(ANN$maxit )) { opt@ANN$maxit <- ANN$maxit }
  }

  if(!is.null(SRE)){
    if(!is.null(SRE$quant )) { opt@SRE$quant <- SRE$quant }
  }

  if(!is.null(FDA)){
#     if(!is.null(FDA$type )) { opt@FDA$type <- FDA$type }
#     if(!is.null(FDA$interaction.level )) { opt@FDA$interaction.level <- FDA$interaction.level }
    if(!is.null(FDA$method )) { opt@FDA$method <- FDA$method }  
  }

  if(!is.null(MARS)){
    if(!is.null(MARS$degree )) { opt@MARS$degree <- MARS$degree }
    if(!is.null(MARS$penalty )) { opt@MARS$penalty <- MARS$penalty }
    if(!is.null(MARS$thresh )) { opt@MARS$thresh <- MARS$thresh }
    if(!is.null(MARS$prune )) { opt@MARS$prune <- MARS$prune }
  }

  if(!is.null(RF)){
    if(!is.null(RF$type )) { opt@RF$type <- RF$type }
#     if(!is.null(RF$interaction.level )) { opt@RF$interaction.level <- RF$interaction.level }
#     if(!is.null(RF$do.classif )) { opt@RF$do.classif <- RF$do.classif }
    if(!is.null(RF$ntree )) { opt@RF$ntree <- RF$ntree }
    if(!is.null(RF$mtry )) { opt@RF$mtry <- RF$mtry }
    if(!is.null(RF$nodesize )) { opt@RF$nodesize <- RF$nodesize }
    if(!is.null(RF$maxnodes )) { opt@RF$maxnodes <- RF$maxnodes }
  }

  if(!is.null(MAXENT)){
    if(!is.null(MAXENT$path_to_maxent.jar )) {
      opt@MAXENT$path_to_maxent.jar <- normalizePath(sub("maxent.jar", "", MAXENT$path_to_maxent.jar)) # ensure path format validity
      } else {opt@MAXENT$path_to_maxent.jar <- getwd()}
    if(!is.null(MAXENT$maximumiterations )) { opt@MAXENT$maximumiterations <- MAXENT$maximumiterations }
    if(!is.null(MAXENT$visible )) { opt@MAXENT$visible <- MAXENT$visible }
    if(!is.null(MAXENT$linear )) { opt@MAXENT$linear <- MAXENT$linear }
    if(!is.null(MAXENT$quadratic )) { opt@MAXENT$quadratic <- MAXENT$quadratic }
    if(!is.null(MAXENT$product )) { opt@MAXENT$product <- MAXENT$product }
    if(!is.null(MAXENT$threshold )) { opt@MAXENT$threshold <- MAXENT$threshold }
    if(!is.null(MAXENT$hinge )) { opt@MAXENT$hinge <- MAXENT$hinge }
    if(!is.null(MAXENT$lq2lqptthreshold )) { opt@MAXENT$lq2lqptthreshold <- MAXENT$lq2lqptthreshold }
    if(!is.null(MAXENT$l2lqthreshold )) { opt@MAXENT$l2lqthreshold <- MAXENT$l2lqthreshold }
    if(!is.null(MAXENT$hingethreshold )) { opt@MAXENT$hingethreshold <- MAXENT$hingethreshold }
    if(!is.null(MAXENT$beta_threshold )) { opt@MAXENT$beta_threshold <- MAXENT$beta_threshold }
    if(!is.null(MAXENT$beta_categorical )) { opt@MAXENT$beta_categorical <- MAXENT$beta_categorical }
    if(!is.null(MAXENT$beta_lqp )) { opt@MAXENT$beta_lqp <- MAXENT$beta_lqp }
    if(!is.null(MAXENT$beta_hinge )) { opt@MAXENT$beta_hinge <- MAXENT$beta_hinge }
    if(!is.null(MAXENT$defaultprevalence )) { opt@MAXENT$defaultprevalence <- MAXENT$defaultprevalence }
  } else{
    opt@MAXENT$path_to_maxent.jar <- getwd()
  }
  
  test <- as.logical(validObject(object = opt, test = TRUE, complete = FALSE))
  
  if(!test){
    cat("\n\n!!! NULL object returned because of invalid parameters given !!!")
    return(NULL)
  }

  return(opt)
}

Print_Default_ModelingOptions <- function(){
  cat('\n Defaut modeling options. copy, change what you want paste it as arg to BIOMOD_ModelingOptions\n\n')
  
  opt_tmp <- BIOMOD_ModelingOptions()
  print(opt_tmp)
}
