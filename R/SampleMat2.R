SampleMat2 <-
function(ref, ratio)
{
    ntot <- length(ref)
    npres<- sum(ref)    
    ncal <- ceiling(ntot*ratio)

    pres <- sample(which(ref==1), ceiling(npres*ratio))
    absc <- sample(which(ref==0), ncal-length(pres))
    
    samprows <- list("calibration"=c(pres,absc), "evaluation"=(1:ntot)[-c(pres,absc)])
    return(samprows)
}

