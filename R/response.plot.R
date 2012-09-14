`response.plot` <-
function(model, Data, show.variables=seq(1:ncol(Data)), save.file="no", name="response_curve", ImageSize=480, plot=TRUE){
  
  cat("\n! Deprecated function, please use response.plot2 instead !")
  return(TRUE)
#   if(inherits(Data,"Raster")){
#     cat("\n   > Extracting raster infos..")
#     DataTmp <- matrix(0,ncol=nlayers(Data), nrow=100)
#     colnames(DataTmp) <- names(Data)
#     maxVal <- maxValue(Data)
#     minVal <- minValue(Data)
#     for(i in 1:ncol(DataTmp)){
#       DataTmp[,i] <- seq(minVal[i],maxVal[i],length.out=100)
#     }
#     Data <- DataTmp
#     rm(list=c('maxVal','minVal','DataTmp'))
#     
#   }
# 
#     if(sum(show.variables > ncol(Data)) > 0) stop("columns wanted in show.variables do not match the data \n")
# 
#     NbVar <- ncol(Data)
#     NbVarShow <- length(show.variables)
#     
#     if(plot==F) temp <- array(0, dim=c(nrow(Data), 2, NbVarShow), dimnames=list(NULL, c("Var", "Pred"), colnames(Data)[show.variables]))
#     
#     
#     #consider if factorial variables :     
#     Xp  <- as.data.frame(matrix(NA, ncol=NbVar, nrow=nrow(Data), dimnames=list(NULL, colnames(Data))))
#     for(i in 1:NbVar){
#         if(is.numeric(Data[,i])) { Xp[,i] <- mean(Data[,i])
#         } 
#         else { 
#         	Xp[, i] <- as.factor(rep(names(which.max(summary(Data[, i]))), nrow(Data)))
#         	levels(Xp[,i]) <- levels(Data[, i])	 	
#         }
#     }   
#       
#     if(substr(class(model)[1],1,4)=="nnet" ) if(sum(search()=="package:nnet")==0) library(nnet)
#     if(class(model)[1]=="rpart") if(sum(search()=="package:rpart")==0) library(rpart)
#     if(class(model)[1]=="mars" | class(model)[1]=="fda") if(sum(search()=="package:mda")==0) library(mda)
#     if("randomForest" %in% class(model)) if(sum(search()=="package:randomForest")==0) library(randomForest,  verbose=FALSE)
#     
#     if(plot){
#     if(save.file=="pdf") pdf(paste(name, "pdf", sep="."))
#     if(save.file=="jpeg") jpeg(paste(name, "jpeg", sep="."), width=ImageSize, height=ImageSize)
#     if(save.file=="tiff") tiff(paste(name, "tiff", sep="."), width=ImageSize, height=ImageSize)
#     if(save.file=="postscript") postscript(paste(name, "eps", sep="."))
#     
#     #plotting window
#     W.width <- ceiling(sqrt(NbVarShow))
#     W.height <- ceiling(NbVarShow/W.width)
#     mat <- matrix(c(rep(1,W.width), 1:(W.height*W.width)+1), ncol=W.width, byrow=TRUE) 
#     layout(mat, widths=rep(1,W.width), heights=c(0.3,rep(1,W.height)))
#     
#     par(mar = c(0.1, 0.1, 0.1, 0.1))
#     plot(x=c(-1,1),y=c(0,1),xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE)
#     polygon(x=c(-2,-2,2,2),y=c(-2,2,2,-2),col="#f5fcba",border=NA)
#     text(x=0.5, y=0.8, pos=1, cex=1.6, labels=paste("Response curves ", class(model)[1], sep=""),col="#4c57eb")
#     par(mar = c(2,2,3.5,1))
# 
# 	}
#     for(i in 1:NbVar){ if(sum(i==show.variables) > 0){
#     
#             #consider if factorial variables :
#             if(!is.factor(Data[,i])){  
#                 xr <- range(Data[,i])
#                 Xp1 <- Xp
#                 Xp1[,i] <- seq(xr[1], xr[2],  len=nrow(Data))
#             } else {
#                 Xp1 <- Xp
#                 Nrepcat <- floor(nrow(Data)/length(levels(Data[,i])))
#                 Xp1[,i] <- as.factor(c(rep(levels(Data[,i])[1], nrow(Data)-(Nrepcat*length(levels(Data[,i])))), rep(levels(Data[,i]), each=Nrepcat)))
#         
#             }
#         
#             if(class(model)[1]=="glm" | class(model)[1]=="gam") Xf <- predict(model, as.data.frame(Xp1), type="response")
#             if(class(model)[1]=="gbm") Xf <-  predict.gbm(model, as.data.frame(Xp1), model$n.trees, type="response") 
#             if(class(model)[1]=="rpart") Xf <- as.numeric(predict(model, Xp1, type="vector"))
#             if(substr(class(model)[1],1,4)=="nnet" ) Xf <- as.numeric(predict(model, as.data.frame(Xp1), type="raw"))
#             if(class(model)[1]=="mars") Xf <- as.numeric(predict(model, as.data.frame(Xp1)))
#             if(class(model)[1]=="fda") Xf <- predict(model, as.data.frame(Xp1), type="post")[,2]
#             if("randomForest" %in% class(model)) Xf <- predict(model, as.data.frame(Xp1), type="prob")[,2]
#       
#       
#             #rescaling preds (not possible to use rescaling_GLM -> no info on calib data)
#             if(class(model)[1]=="mars" | substr(class(model)[1],1,4)=="nnet"  | class(model)[1]=="fda" ){ 
#                 OriMinMax <- range(Xf)	
#                 Xf <- (Xf - min(OriMinMax)) / (max(OriMinMax)-min(OriMinMax))
#                 Xf[Xf<0]<-0
#                 Xf[Xf>1]<-1            
#             }
# #             cat("no")
# 			if(plot) {
# 				plot(Xp1[ ,i], Xf, ylim=c(0,1), xlab="", ylab="", type="l", main=names(Data)[i])
# 				rug(Data[ ,i])
# 			}	
# 			else{ 
# 				temp[,1,i] <-Xp1[ ,i]; temp[,2,i] <- Xf
# 			}     
#     }}# i loop for variables
#    
#    
#     if(save.file=="pdf" | save.file=="jpeg" | save.file=="tiff" | save.file=="postscript") dev.off()
#     if(plot==F) return(temp)
#    
#   #  if(substr(class(model)[1],1,4)=="nnet" )  detach(package:nnet)
#    # if(class(model)[1]=="rpart") detach(package:rpart)
#    # if(class(model)[1]=="mars" | class(model)[1]=="fda") detach(package:mda)
#    # if(class(model)[1]=="randomForest") detach(package:randomForest)            
} 


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
`response.plot2` <- function(models, 
                             Data, 
                             show.variables=seq(1:ncol(Data)),
                             do.bivariate = FALSE,
                             fixed.var.metric = 'mean',
                             save.file="no", 
                             name="response_curve", 
                             ImageSize=480, 
                             plot=TRUE,
                             ...){
  
  # 1. args checking -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  add.args <- list(...)
  args <- .response.plot2.check.arg(models, Data, show.variables, save.file, name, ImageSize, plot, fixed.var.metric, do.bivariate, add.args)
  
  models <- args$models
  Data <- args$Data
  show.variables <- args$show.variables
  save.file <- args$save.file
  name <- args$name
  ImageSize <- args$ImageSize
  plot <- args$plot
  fixed.var.metric <- args$fixed.var.metric
  do.bivariate <- args$do.bivariate
  nb.pts <- args$nb.pts

  # 2. build function outputs
  
  ## array for monavariate res
  array.mono.out <- array(0, 
                     dim=c(nb.pts,2,length(show.variables), length(models)), 
                     dimnames=list(NULL,  c("Var", "Pred"), show.variables, models) )
  if(do.bivariate){
    ## array for bivariate res
    dim3.array.bi.out <- c()
    for(i in 1:(length(show.variables)-1)){
      for(j in (i+1):length(show.variables)){
        dim3.array.bi.out <- c(dim3.array.bi.out, paste(show.variables[i],show.variables[j],sep="-"))
      }
    }
    
    array.bi.out <- array(0, 
                       dim=c(nb.pts,3,length(dim3.array.bi.out), length(models)), 
                       dimnames=list(NULL,  c("Var1", "Var2", "Pred"), dim3.array.bi.out, models) )    
  }

  # Create a ranged data table
  Data.r <- data.frame(matrix(NA, nrow=nb.pts, ncol=ncol(Data)))
  colnames(Data.r) <- colnames(Data)
  for(i in 1:ncol(Data)){
    if(is.numeric(Data[,i])){
      Data.r[,i] <- switch(fixed.var.metric,
                           mean = rep(mean(Data[,i]), nb.pts),
                           median = rep(median(Data[,i]), nb.pts),
                           min = rep(min(Data[,i]), nb.pts),
                           max = rep(max(Data[,i]), nb.pts))
    } else{
      Data.r[,i] <- switch(fixed.var.metric,
                           mean = rep(levels(as.factor(Data[,i]))[round(length(levels(as.factor(Data[,i]))) / 2 )], nb.pts),
                           median = rep(levels(as.factor(Data[,i]))[round(length(levels(as.factor(Data[,i]))) / 2 )], nb.pts),
                           min = rep(levels(as.factor(Data[,i]))[1], nb.pts),
                           max = rep(levels(as.factor(Data[,i]))[length(levels(as.factor(Data[,i])))], nb.pts))
      Data.r[,i] <- as.factor(Data.r[,i])
    }
  }
  
  if(plot){
    # X. Open a graphic file for plotting restults
    if(save.file=="pdf") pdf(paste(name, "pdf", sep="."))
    if(save.file=="jpeg") jpeg(paste(name, "jpeg", sep="."), width=ImageSize, height=ImageSize)
    if(save.file=="tiff") tiff(paste(name, "tiff", sep="."), width=ImageSize, height=ImageSize)
    if(save.file=="postscript") postscript(paste(name, "eps", sep="."))
    
    # XX. parametrize our plot window
      
    if(!do.bivariate){
      nb.graphs <- length(show.variables)
    } else{
      nb.graphs <- length(models) *  ( (length(show.variables)-1) * length(show.variables) / 2 )
    }
    
    W.width <- ceiling(sqrt(nb.graphs))
    W.height <- ceiling(nb.graphs/W.width)    
    
    mat <- matrix(c(rep(1,W.width), 1:(W.height*W.width)+1), ncol=W.width, byrow=TRUE) 
    layout(mat, widths=rep(1,W.width), heights=c(0.3,rep(1,W.height)))
    
    par(mar = c(0.1, 0.1, 0.1, 0.1))
    plot(x=c(-1,1),y=c(0,1),xlim=c(0,1),ylim=c(0,1),type="n",axes=FALSE)
    polygon(x=c(-2,-2,2,2),y=c(-2,2,2,-2),col="#f5fcba",border=NA)
    text(x=0.5, y=0.8, pos=1, cex=1.6, labels=paste("Response curves for ", .extractModelNamesInfo(models[1],"species"), "'s ", .extractModelNamesInfo(models[1],"models"),sep=""),  ,col="#4c57eb")
    par(mar = c(2,2,3.5,1))      
  } 

  
  
  if(!do.bivariate){  
    for(vari in show.variables){
    	if(plot) {
  #       frame()
        plot(0,0,col="white",xlim=c(min(Data[,vari]), max(Data[,vari])), ylim=c(0,1), main=vari, ann=TRUE, bty="o",xaxs="r", xaxt="s")
  			rug(Data[ ,vari])
  		}
      
      for(model in models){
        
        # 0. get model
        mod <- get(model)
        
        # 1. load library if some needed
        if(substr(class(mod)[1],1,4)=="nnet" ) if(sum(search()=="package:nnet")==0) library(nnet)
        if(class(mod)[1]=="rpart") if(sum(search()=="package:rpart")==0) library(rpart)
        if(class(mod)[1]=="mars" | class(mod)[1]=="fda") if(sum(search()=="package:mda")==0) library(mda)
        if("randomForest" %in% class(mod)) if(sum(search()=="package:randomForest")==0) library(randomForest,  verbose=FALSE)
        if(inherits(mod, 'gbm')) if(sum(search()=="package:gbm")==0) library(gbm,  verbose=FALSE)
          
        # 2. do projections
        pts.tmp <- seq(min(Data[,vari]), max(Data[,vari]), length.out=nb.pts)
        
        Data.r.tmp <- Data.r
        Data.r.tmp[,vari] <- pts.tmp
        
        if(inherits(mod,'nnet')){ set.seed(555); proj.tmp <- as.numeric(predict(mod, Data.r.tmp, type = "raw")) }
        if(inherits(mod,'rpart')){ proj.tmp <- as.numeric(predict(mod, Data.r.tmp, type="prob")[,2]) }
        if(inherits(mod,'gam')){ proj.tmp <- predict(object=mod, newdata=Data.r.tmp, type="response") }
        if(inherits(mod,'gbm')){ best.iter <- gbm.perf(mod, method = "cv", plot.it = FALSE); proj.tmp <- predict.gbm(mod, Data.r.tmp, best.iter, type = "response") }
        if(inherits(mod,'glm') & !inherits(mod,'gam')){ proj.tmp <- .testnull(mod, Prev=0.5,  Data.r.tmp) }
        if(inherits(mod,'fda')){ proj.tmp <- as.numeric(predict(mod, Data.r.tmp, type = "posterior")[, 2]) }
        if(inherits(mod,'mars')){ proj.tmp <- as.numeric(predict(mod, Data.r.tmp)) }
        if(inherits(mod,'randomForest')){ proj.tmp <- as.numeric(predict(mod,Data.r.tmp, type='prob')[,'1']) }
        if(inherits(mod,'maxent_model')){ proj.tmp <- as.numeric(predict(object=mod, newdata=Data.r.tmp, proj_name='RespPlotTmp', rm_tmp_files=TRUE)) }
        
        # 3. Rescaling stuff
        ## TO DO
        
        # 4. Ploting results
    		if(plot) {
  				lines(pts.tmp, proj.tmp)
  			}
        
        # 5. Storing results
        array.mono.out[,"Var",vari,model] <- pts.tmp
        array.mono.out[,"Pred",vari,model] <- proj.tmp
      }    
      
    }
  } else{ ## bivariate case
    for(vari1 in show.variables[-length(show.variables)]){
      for(vari2 in show.variables[-(1:which(show.variables == vari1))]){
#         if(plot) {
#     #       frame()
#           persp(0,0,col="white",xlim=c(min(Data[,vari]), max(Data[,vari])), ylim=c(0,1), main=vari, ann=TRUE, bty="o",xaxs="r", xaxt="s")
#     			rug(Data[ ,vari])
#     		}
        
        for(model in models){
          
          # 0. get model
          mod <- get(model)
          
          # 1. load library if some needed
          if(substr(class(mod)[1],1,4)=="nnet" ) if(sum(search()=="package:nnet")==0) library(nnet)
          if(class(mod)[1]=="rpart") if(sum(search()=="package:rpart")==0) library(rpart)
          if(class(mod)[1]=="mars" | class(mod)[1]=="fda") if(sum(search()=="package:mda")==0) library(mda)
          if("randomForest" %in% class(mod)) if(sum(search()=="package:randomForest")==0) library(randomForest,  verbose=FALSE)
          if(inherits(mod, 'gbm')) if(sum(search()=="package:gbm")==0) library(gbm,  verbose=FALSE)
            
          # 2. do projections
          pts.tmp1 <- rep(seq(min(Data[,vari1]), max(Data[,vari1]), length.out=sqrt(nb.pts)),each=sqrt(nb.pts))
          pts.tmp2 <- rep(seq(min(Data[,vari2]), max(Data[,vari2]), length.out=sqrt(nb.pts)),sqrt(nb.pts))
          
          Data.r.tmp <- Data.r
          Data.r.tmp[,vari1] <- pts.tmp1
          Data.r.tmp[,vari2] <- pts.tmp2
          
          if(inherits(mod,'nnet')){ set.seed(555); proj.tmp <- as.numeric(predict(mod, Data.r.tmp, type = "raw")) }
          if(inherits(mod,'rpart')){ proj.tmp <- as.numeric(predict(mod, Data.r.tmp, type="prob")[,2]) }
          if(inherits(mod,'gam')){ proj.tmp <- predict(object=mod, newdata=Data.r.tmp, type="response") }
          if(inherits(mod,'gbm')){ best.iter <- gbm.perf(mod, method = "cv", plot.it = FALSE); proj.tmp <- predict.gbm(mod, Data.r.tmp, best.iter, type = "response") }
          if(inherits(mod,'glm') & !inherits(mod,'gam')){ proj.tmp <- .testnull(mod, Prev=0.5,  Data.r.tmp) }
          if(inherits(mod,'fda')){ proj.tmp <- as.numeric(predict(mod, Data.r.tmp, type = "posterior")[, 2]) }
          if(inherits(mod,'mars')){ proj.tmp <- as.numeric(predict(mod, Data.r.tmp)) }
          if(inherits(mod,'randomForest')){ proj.tmp <- as.numeric(predict(mod,Data.r.tmp, type='prob')[,'1']) }
          if(inherits(mod,'maxent_model')){ proj.tmp <- as.numeric(predict(object=mod,newdata=Data.r.tmp,xy=data.frame(x=rep(0,nrow(Data.r.tmp)), y=rep(0,nrow(Data.r.tmp))), proj_name='RespPlotTmp', rm_tmp_files=TRUE)) }
          
          # 3. Rescaling stuff
          ## TO DO
          
          # 4. Storing results
          array.bi.out[,"Var1",paste(vari1,vari2,sep="-"),model] <- pts.tmp1
          array.bi.out[,"Var2",paste(vari1,vari2,sep="-"),model] <- pts.tmp2
          array.bi.out[,"Pred",paste(vari1,vari2,sep="-"),model] <- proj.tmp
          
          # 5. Ploting results
        	if(plot) {
            # reformating results to perform a persp plot
            pts.tmp1 <- sort(unique(pts.tmp1))
            pts.tmp2 <- sort(unique(pts.tmp2))
            proj.tmp <- matrix(proj.tmp, ncol=length(pts.tmp2), byrow=FALSE)
            
            # build color scale
            ncz <- length(pts.tmp2)
            nrz <- length(pts.tmp1)
            # Create a function interpolating colors in the range of specified colors
            jet.colors <- colorRampPalette(c("red", "orange", "green"))
            # Generate the desired number of colors from this palette
            nbcol <- 50
            color <- jet.colors(nbcol)
            # Compute the z-value at the facet centres
            zfacet <- proj.tmp[-1, -1] + proj.tmp[-1, -ncz] + proj.tmp[-nrz, -1] + proj.tmp[-nrz, -ncz]
            # Recode facet z-values into color indices
            facetcol <- cut(zfacet, nbcol)
            
    				persp(x=pts.tmp1,y=pts.tmp2,z=proj.tmp, xlab = vari1, ylab=vari2, zlab="pred", theta = 30, phi = 30,
              expand = 0.5, col = color[facetcol], ltheta = 120, shade = 0.25, ticktype = "simple", main = paste(.extractModelNamesInfo(model,"data.set"), " ", .extractModelNamesInfo(model,"run.eval") ,sep=""), cex.main = 0.9, cex.axis=0.7)
    			}
          
        }        
      }
    }
  }

  # XXX. Close file
  if(save.file=="pdf" | save.file=="jpeg" | save.file=="tiff" | save.file=="postscript") dev.off()
  
  # delete temp files if somes has been created
  if(file.exists(file.path(.extractModelNamesInfo(models[[1]],"species"),'RespPlotTmp'))){
     unlink(path.expand(file.path(.extractModelNamesInfo(models[[1]],"species"),'RespPlotTmp')), recursive=TRUE, force=TRUE)
  }
  
  if(!do.bivariate){
    invisible(array.mono.out)
  } else{
    invisible(array.bi.out)
  }
    
}
 
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
.response.plot2.check.arg <- function(models, Data, show.variables, save.file, name, ImageSize, plot, fixed.var.metric, do.bivariate, add.args){
  
  # 1. check add.args
  if(sum(! (names(add.args) %in% c("nb.pts","xy"))) > 0){
    warning(paste(toString(names(add.args)[which(! (names(add.args) %in% c("nb.pts")))]), " are unknown arguments", sep="" ))
  }
  
  
  ### check of models args =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  if(!is.character(models)){
    stop("models must be a character vector of models names")
  }
  for(mod in models){
    if(!exists(mod)){
      stop("you need to load the models selected!")
    }
  }
  
  ### defining the number split in each variables range =-=-=-=-=- #
  if(!is.null(add.args$nb.pts)){
    if(do.bivariate){
      # total number of points is the square of the difined
      add.args$nb.pts <- add.args$nb.pts^2
    }
  } else{
    if(!do.bivariate){
      add.args$nb.pts <- 100
    } else{
      add.args$nb.pts <- 25^2
    }
  }
  
  ### check of data args =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  if(inherits(Data,"Raster")){
    cat("\n   > Extracting raster infos..")
    DataTmp <- matrix(0,ncol=nlayers(Data), nrow=add.args$nb.pts)
    colnames(DataTmp) <- names(Data)
    maxVal <- maxValue(Data)
    minVal <- minValue(Data)
    for(i in 1:ncol(DataTmp)){
      DataTmp[,i] <- seq(minVal[i],maxVal[i],length.out=add.args$nb.pts)
    }
    Data <- DataTmp
    rm(list=c('maxVal','minVal','DataTmp'))
    
  }
    
  ### check show.variables arg -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if( ( length(show.variables) > ncol(Data) ) | (sum(!(show.variables %in% colnames(Data)))) ) stop("columns wanted in show.variables do not match the data \n")
    
  ### check save.file arg -=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

  


  
  
  # TO DO 
  return(list(models = models,
              Data = Data,
              show.variables = show.variables, 
              save.file = save.file, 
              name = name, 
              ImageSize = ImageSize, 
              plot = plot,
              fixed.var.metric = fixed.var.metric,
              do.bivariate = do.bivariate,
              nb.pts = add.args$nb.pts))
}

###
# list of supported additional arguments  
# - nb.pts  

