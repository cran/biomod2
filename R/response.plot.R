`response.plot` <-
function(model, Data, show.variables=seq(1:ncol(Data)), save.file="no", name="response_curve", ImageSize=480, plot=TRUE){
  
  cat("\n! Deprecated function, please use response.plot2 instead !")
  return(TRUE)      
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
  on_0_1000 <- add.args$on_0_1000
  
  if(is.null(on_0_1000)) on_0_1000 <- FALSE
  
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
        if(on_0_1000) ylim <- c(0,1000) else ylim <- c(0,1) 

        plot(0,0,col="white",xlim=c(min(Data[,vari]), max(Data[,vari])), ylim=ylim, main=vari, ann=TRUE, bty="o",xaxs="r", xaxt="s")
  			rug(Data[ ,vari])
  		}
      
      for(model in models){
        
        # 0. get model
        mod <- get(model)

        # 2. do projections
        pts.tmp <- seq(min(Data[,vari]), max(Data[,vari]), length.out=nb.pts)
        
        Data.r.tmp <- Data.r
        Data.r.tmp[,vari] <- pts.tmp
        
        proj.tmp <- predict(mod, Data.r.tmp, on_0_1000=on_0_1000)

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
      
        for(model in models){
          
          # 0. get model
          mod <- get(model)
            
          # 2. do projections
          pts.tmp1 <- rep(seq(min(Data[,vari1]), max(Data[,vari1]), length.out=sqrt(nb.pts)),each=sqrt(nb.pts))
          pts.tmp2 <- rep(seq(min(Data[,vari2]), max(Data[,vari2]), length.out=sqrt(nb.pts)),sqrt(nb.pts))
          
          Data.r.tmp <- Data.r
          Data.r.tmp[,vari1] <- pts.tmp1
          Data.r.tmp[,vari2] <- pts.tmp2
          
          proj.tmp <- predict(mod, Data.r.tmp, on_0_1000=on_0_1000)
          
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

