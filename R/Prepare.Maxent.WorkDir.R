.Prepare.Maxent.WorkDir <- function(Data, xy, calibLines=NULL, RunName=NULL, VarImport=0, evalData=NULL, evalxy=NULL, species.name=NULL, modeling.id='' ){
  cat('\n\tCreating Maxent Temp Proj Data..')
  
  ## initialise output 
  MWD <- list()
  class(MWD) <- "maxent_workdir_info"
  
  ## default parameters setting
  if(is.null(RunName)) RunName <- colnames(Data)[1]
  if(is.null(species.name)) species.name <- colnames(Data)[1]
  if(is.null(calibLines)) calibLines <- rep(T,nrow(Data))
  
  ## define all paths to files needed by MAXENT
  m_workdir <- file.path(species.name,'models',modeling.id,paste('m_',sub(".","",as.character(format(Sys.time(), "%OS6")), fixed=T),sep=""))
  # check wordir unicity
  while(file.exists(m_workdir)){
    m_workdir <- file.path(species.name,'models',modeling.id,paste('m_',sub(".","",as.character(format(Sys.time(), "%OS6")), fixed=T),sep=""))
  }
  MWD$m_workdir <- m_workdir
  
  m_outdir <- file.path(species.name,'models',modeling.id, paste(RunName,'_MAXENT_outputs',sep=''))
  MWD$m_outdir <- m_outdir
  
  ## directories creation
  dir.create(m_workdir, showWarnings=FALSE, recursive=TRUE, mode='777')
  dir.create(m_outdir, showWarnings=FALSE, recursive=TRUE, mode='777')
  
  # Presences Data
  m_speciesFile <- file.path(m_workdir,"Sp_swd.csv")
  MWD$m_speciesFile <- m_speciesFile
  presLines <- which((Data[,1]==1) & calibLines)
  absLines <- which((Data[,1]==0) & calibLines)
  Sp_swd <- cbind(rep(RunName,length(presLines)),
                      xy[presLines,],
                      Data[presLines,2:ncol(Data),drop=FALSE])
  colnames(Sp_swd) <- c('species','X','Y',colnames(Data)[2:ncol(Data)])
  write.table(Sp_swd, file=m_speciesFile,  quote=FALSE, row.names=FALSE, sep=",")
  
  # Background Data
  m_backgroundFile <- file.path(m_workdir,"Back_swd.csv")
  MWD$m_backgroundFile <- m_backgroundFile
  # keep only 0 of calib lines
  Back_swd <- cbind(rep("background",length(absLines)),xy[absLines,],Data[absLines,2:ncol(Data),drop=FALSE])
  colnames(Back_swd)  <- c("background",colnames(Back_swd)[-1])
  write.table(Back_swd, file=m_backgroundFile, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
  
  # Prediction Data
  m_predictDir <- file.path(m_workdir,"Predictions")
  dir.create(m_predictDir, showWarnings=FALSE, recursive=TRUE)
  
  m_predictFile <- file.path(m_predictDir, "Pred_swd.csv")
  MWD$m_predictFile <- m_predictFile
  Pred_swd <- cbind(rep("predict",nrow(xy)),xy,Data[,2:ncol(Data),drop=FALSE])
  colnames(Pred_swd)  <- c("predict",colnames(Back_swd)[-1])
  write.table(Pred_swd, file=m_predictFile, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
  
  ### USELESS ###
  # dealing with variable importances stuff
#   if( VarImport > 0){
#     for( vari in colnames(Data)[-1] )
#       for(vi in 1:VarImport){
#         proj_tmp <- Pred_swd
#         proj_tmp[,1] <- rep(paste(vari,'_',vi,sep=""),nrow(proj_tmp))
#         proj_tmp[,vari] <- sample(proj_tmp[,vari])
#         write.table(proj_tmp, file=file.path(species.name,'models',modeling.id,"MaxentTmpData","Pred",paste(vari,'_',vi,"_swd.csv",sep="")), quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
#       }
#   }
  
  # dealing with independent evaluation data
  if(!is.null(evalData)){
    m_predictEvalFile <- file.path(m_predictDir, "PredEval_swd.csv")
    MWD$m_predictEvalFile <- m_predictEvalFile
    Pred_eval_swd <- cbind(rep("predictEval",nrow(evalxy)),evalxy,evalData[,2:ncol(evalData),drop=FALSE])
    colnames(Pred_eval_swd)  <- c("predict",colnames(Back_swd)[-1])
    write.table(Pred_eval_swd, file=m_predictEvalFile, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
  }
  
  return(MWD)
}

.Delete.Maxent.WorkDir <- function(MWD, silent=FALSE){
  if(!silent) cat('\n\tRemoving Maxent Temp Data..')
  if(class(MWD) == "maxent_workdir_info"){
    unlink(MWD$m_workdir,recursive = TRUE, force = TRUE)
  } else{
    if(!silent) cat('\n\t! Invalid maxent work dir object -> MAXENT temp files have not been removed')
  }
}

# Maxent Projection working directory preparation -=-=-=-=-=-=-=- #

setGeneric(".Prepare.Maxent.Proj.WorkDir", 
            def = function(Data, ...){
              standardGeneric( ".Prepare.Maxent.Proj.WorkDir" )
            } )


setMethod('.Prepare.Maxent.Proj.WorkDir', signature(Data='data.frame'),
          def = function(Data, xy, species.name =".", proj.name=".", silent=FALSE){
            ## initialise output 
            MWD <- list()
            class(MWD) <- "maxent_workdir_info"
            
            if(!silent) cat('\n\t\tCreating Maxent Temp Proj Data...')
            if(is.null(xy)) xy <- matrix(1,nrow=nrow(Data), ncol=2, dimnames=list(NULL, c("X","Y")))

            ## define all paths to files needed by MAXENT
            m_workdir <- file.path(species.name,proj.name,paste('m_',sub(".","",as.character(format(Sys.time(), "%OS6")), fixed=T),sep=""))
            # check wordir unicity
            while(file.exists(m_workdir)){
              m_workdir <- file.path(species.name,proj.name,paste('m_',sub(".","",as.character(format(Sys.time(), "%OS6")), fixed=T),sep=""))
            }
            MWD$m_workdir <- m_workdir
                      
            dir.create(m_workdir, recursive=TRUE, showWarnings=FALSE)
            
            # Proj Data
            m_predictFile <- file.path(m_workdir, "Pred_swd.csv")
            MWD$m_predictFile <- m_predictFile
            Proj_swd <- cbind(rep("proj",nrow(xy)),xy,Data)
            colnames(Proj_swd)  <- c("proj","X","Y",colnames(Data))
            write.table(Proj_swd, file=m_predictFile, quote=FALSE, row.names=FALSE, col.names=TRUE, sep=",")
            
            return(MWD)
            })


setMethod('.Prepare.Maxent.Proj.WorkDir', signature(Data='RasterStack'),
          def = function(Data, species.name =".",proj.name=".", silent=FALSE){
            ## initialise output 
            MWD <- list()
            class(MWD) <- "maxent_workdir_info"
            
            if(!silent) cat('\n\t\tCreating Maxent Temp Proj Data...')
            
            ## define all paths to files needed by MAXENT
            m_workdir <- file.path(species.name,proj.name,paste('m_',sub(".","",as.character(format(Sys.time(), "%OS6")), fixed=T),sep=""))
            # check wordir unicity
            while(file.exists(m_workdir)){
              m_workdir <- file.path(species.name,proj.name,paste('m_',sub(".","",as.character(format(Sys.time(), "%OS6")), fixed=T),sep=""))
            }
            MWD$m_workdir <- m_workdir
            dir.create(m_workdir, recursive=TRUE, showWarnings=FALSE)
            
            # Proj Data
            m_predictFile <- NULL
            for(l in names(Data)){
              m_predictFileTmp <- file.path(m_workdir,paste(l,'.asc',sep='')) 
              if(! file.exists(m_predictFileTmp)){
                if(!silent) cat("\n\t\t\t>",l ,"\t:\t" )
                if(grepl(".asc", filename(raster::subset(Data,l,drop=TRUE)) ) ){
                  if(!silent) cat("coping ascii file")
                  file.copy(filename(raster::subset(Data,l,drop=TRUE)), m_predictFileTmp)
                } else{
                  if(!silent) cat("creating ascii file")
                  writeRaster(raster::subset(Data,l,drop=TRUE), filename=m_predictFileTmp,
                              format='ascii', overwrite=TRUE)        
                }
                
              } else{
                if(!silent) cat("\n", m_predictFileTmp ,'already created !')
              }
              m_predictFile <- c(m_predictFile, m_predictFileTmp)
            }
            MWD$m_predictFile <- m_predictFile
            
            return(MWD)
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

#       if(grepl(".asc", filename(raster::subset(Data,l,drop=TRUE)) ) ){
#         cat("\n coping ascii file")
#         file.copy(filename(raster::subset(Data,l,drop=TRUE)), file.path(proj_name,'MaxentTmpData', 'Proj' ,paste(l,'.asc',sep='')))
#       } else{
#         cat("\n creating ascii file")
#         writeRaster(raster::subset(Data,l,drop=TRUE), filename=file.path(proj_name,'MaxentTmpData', 'Proj' ,paste(l,'.asc',sep='')),
#             format='ascii', overwrite=TRUE)        
#       }
# 
#     } else{
#       cat("\n", file.path(proj_name,'MaxentTmpData',paste(l,'.asc',sep='')),'already created !')
#     }
# 
#   }
# }
