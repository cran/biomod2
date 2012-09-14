.Prepare.Maxent.WorkDir <- function(Data, xy, calibLines, RunName=NULL, VarImport=0, evalData=NULL, evalxy=NULL, species.name=NULL ){
  cat('\n\tCreating Maxent Temp Proj Data..')
  if(is.null(RunName)) RunName <- colnames(Data)[1]
  if(is.null(species.name)) species.name <- colnames(Data)[1]
  
  dir.create(paste(getwd(),'/',species.name,'/MaxentTmpData', sep=""), showWarnings=FALSE, recursive=TRUE)
  dir.create(paste(getwd(),'/',species.name,'/models/',RunName,'_MAXENT_outputs',sep=''), showWarnings=FALSE, recursive=TRUE)
  
  # Presences Data
  presLines <- which((Data[,1]==1) & calibLines)
  absLines <- which((Data[,1]==0) & calibLines)
  Sp_swd <- cbind(rep(RunName,length(presLines)),
                      xy[presLines,],
                      Data[presLines,2:ncol(Data)])
  colnames(Sp_swd) <- c('specie','X','Y',colnames(Data)[2:ncol(Data)])
  write.table(Sp_swd, file=paste(getwd(),'/',species.name,"/MaxentTmpData/Sp_swd.csv",sep=""), quote=FALSE, row.names=FALSE, sep=",")
  
  # Background Data
  # keep only 0 of calib lines
  Back_swd <- cbind(rep("background",length(absLines)),xy[absLines,],Data[absLines,2:ncol(Data)])
  colnames(Back_swd)  <- c("background",colnames(Back_swd)[-1])
  write.table(Back_swd, file=paste(getwd(),'/',species.name,"/MaxentTmpData/Back_swd.csv",sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
  
  # Prediction Data
  dir.create(paste(getwd(),'/',species.name,'/MaxentTmpData/Pred', sep=""), showWarnings=FALSE)
  
  Pred_swd <- cbind(rep("predict",nrow(xy)),xy,Data[,2:ncol(Data)])
  colnames(Pred_swd)  <- c("predict",colnames(Back_swd)[-1])
  write.table(Pred_swd, file=paste(getwd(),'/',species.name,"/MaxentTmpData/Pred/Pred_swd.csv",sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
  
  # dealing with variable importances stuff
  if( VarImport > 0){
    for( vari in colnames(Data)[-1] )
      for(vi in 1:VarImport){
        proj_tmp <- Pred_swd
        proj_tmp[,1] <- rep(paste(vari,'_',vi,sep=""),nrow(proj_tmp))
        proj_tmp[,vari] <- sample(proj_tmp[,vari])
        write.table(proj_tmp, file=paste(getwd(),'/',species.name,"/MaxentTmpData/Pred/",vari,'_',vi,"_swd.csv",sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
      }
  }
  
  # dealing with independent evaluation data
  if(!is.null(evalData)){
    Pred_eval_swd <- cbind(rep("predictEval",nrow(evalxy)),evalxy,evalData[,2:ncol(evalData)])
    colnames(Pred_eval_swd)  <- c("predict",colnames(Back_swd)[-1])
    write.table(Pred_eval_swd, file=paste(getwd(),'/',species.name,"/MaxentTmpData/Pred/Pred_eval_swd.csv",sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
  }
}

.Delete.Maxent.WorkDir <- function(species.name=NULL){
  cat('\n\tRemoving Maxent Temp Data..')
#   system('rm -rf MaxentTmpData')
  
  if(is.null(species.name)){
    unlink('MaxentTmpData', recursive = TRUE, force = TRUE)
  } else{
     unlink(paste(species.name,"/MaxentTmpData", sep=""), recursive = TRUE, force = TRUE)
  }
  
}

# Maxent Projection working directory preparation -=-=-=-=-=-=-=- #

setGeneric(".Prepare.Maxent.Proj.WorkDir", 
            def = function(Data, ...){
              standardGeneric( ".Prepare.Maxent.Proj.WorkDir" )
            } )


setMethod('.Prepare.Maxent.Proj.WorkDir', signature(Data='data.frame'),
          def = function(Data, xy, proj_name=NULL){
            cat('\n\t\tCreating Maxent Temp Proj Data...')
            
            if(is.null(proj_name)) proj_name <- colnames(Data)[1]
            dir.create(file.path(getwd(),proj_name,'MaxentTmpData'), showWarnings=FALSE, recursive=TRUE)
            
            # Proj Data
            Proj_swd <- cbind(rep("proj",nrow(xy)),xy,Data)
            colnames(Proj_swd)  <- c("proj","X","Y",colnames(Data))
            write.table(Proj_swd, file=file.path(getwd(),proj_name,'MaxentTmpData', 'Proj_swd.csv'), quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
            })


setMethod('.Prepare.Maxent.Proj.WorkDir', signature(Data='RasterStack'),
          def = function(Data, proj_name=NULL){
            cat('\n\t\tCreating Maxent Temp Proj Data...')
            
            if(is.null(proj_name)) proj_name <- colnames(Data)[1]
            dir.create(file.path(getwd(),proj_name,'MaxentTmpData','Proj'), showWarnings=FALSE, recursive=TRUE)
            
            # Proj Data
            for(l in names(Data)){
              if(! file.exists(file.path(proj_name,'MaxentTmpData','Proj',paste(l,'.asc',sep='')))){
                cat("\n\t\t\t>",l ,"\t:\t" )
                if(grepl(".asc", filename(raster:::subset(Data,l,drop=TRUE)) ) ){
                  cat("coping ascii file")
                  file.copy(filename(raster:::subset(Data,l,drop=TRUE)), file.path(proj_name,'MaxentTmpData', 'Proj' ,paste(l,'.asc',sep='')))
                } else{
                  cat("creating ascii file")
                  writeRaster(raster:::subset(Data,l,drop=TRUE), filename=file.path(proj_name,'MaxentTmpData', 'Proj' ,paste(l,'.asc',sep='')),
                              format='ascii', overwrite=TRUE)        
                }
                
              } else{
                cat("\n", file.path(proj_name,'MaxentTmpData','', paste(l,'.asc',sep='')),'already created !')
              }
              
            }
          })

# .Prepare.Maxent.Proj.WorkDir <- function(Data, xy, proj_name=NULL){
#   cat('\n\tCreating Maxent Temp Proj Data..')
#   
#   if(is.null(proj_name)) proj_name <- colnames(Data)[1]
#   dir.create(paste(getwd(),'/',proj_name,'/MaxentTmpData', sep=""), showWarnings=FALSE)
#   
#   # Proj Data
#   Proj_swd <- cbind(rep("proj",nrow(xy)),xy,Data)
#   colnames(Proj_swd)  <- c("proj","X","Y",colnames(Data))
#   write.table(Proj_swd, file=paste(getwd(),'/',proj_name,"/MaxentTmpData/Proj_swd.csv",sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
# }
# 
# .Prepare.Maxent.Proj.Raster.WorkDir <- function(Data, proj_name=NULL){
#   cat('\n\tCreating Maxent Temp Proj Data..')
#   
#   if(is.null(proj_name)){
#     stop("Please refere explicitly a proj name!")
#   }
#   dir.create(paste(getwd(),'/',proj_name,'/MaxentTmpData/Proj', sep=""), showWarnings=FALSE, recursive=TRUE)
#   
#   # Proj Data
#   for(l in names(Data)){
#     if(! file.exists(file.path(proj_name,'MaxentTmpData','Proj',paste(l,'.asc',sep='')))){

#       if(grepl(".asc", filename(raster:::subset(Data,l,drop=TRUE)) ) ){
#         cat("\n coping ascii file")
#         file.copy(filename(raster:::subset(Data,l,drop=TRUE)), file.path(proj_name,'MaxentTmpData', 'Proj' ,paste(l,'.asc',sep='')))
#       } else{
#         cat("\n creating ascii file")
#         writeRaster(raster:::subset(Data,l,drop=TRUE), filename=file.path(proj_name,'MaxentTmpData', 'Proj' ,paste(l,'.asc',sep='')),
#             format='ascii', overwrite=TRUE)        
#       }
# 
#     } else{
#       cat("\n", file.path(proj_name,'MaxentTmpData',paste(l,'.asc',sep='')),'already created !')
#     }
# 
#   }
# }
