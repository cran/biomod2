SampleMat2 <-
function(ref, ratio, as.logi=FALSE)
{
    ntot <- length(ref)
    npres<- sum(ref)    
    ncal <- ceiling(ntot*ratio)

    pres <- sample(which(ref==1), ceiling(npres*ratio))
    absc <- sample(which(ref==0), ncal-length(pres))
    
    if(as.logi){
      calib <- rep(FALSE, ntot)
      calib[c(pres,absc)] <- TRUE
      eval <- !calib
    } else{
      calib <- c(pres,absc)
      eval <- (1:ntot)[-c(pres,absc)]
    }
    
    return(list("calibration"=calib, "evaluation"=eval))
}

