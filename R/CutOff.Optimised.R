`.CutOff.Optimised` <-
function(Obs, Fit){
    
    SumObs <- sum(Obs)
    LengObs <- length(Obs)
    tt <- c(100)
    Cut <- c(0,0,0)
    
    if(length(unique(Fit))==1){
        Cut[1] <- unique(Fit)
        Cut[2] <- 100*sum((Fit>=Cut[1])[Obs==1])/SumObs
        Cut[3] <- 100*sum((Fit<Cut[1])[Obs==0])/(LengObs-SumObs)
        Cut <- t(Cut)
    }
    
    else{
        if(min(Fit)<0) Fit[Fit<0] <- 0
        Quant <- quantile(Fit)
        i <- Quant[1]
        a <- 2
        while(i<=Quant[5]){
            se <- sum((Fit>=i)[Obs==1])/SumObs
            sp <- sum((Fit<i)[Obs==0])/(LengObs-SumObs)
            tt[a] <- abs(se-sp)
            if(tt[a]>tt[a-1]) break
            i <- i+((Quant[5] - Quant[1])/1000)
            a <- a+1
        }
        b <- (i-((Quant[5] - Quant[1])/1000))
        Cut[1] <- b
        Cut[2] <- 100*sum((Fit>=b)[Obs==1])/SumObs
        Cut[3] <- 100*sum((Fit<b)[Obs==0])/(LengObs-SumObs)
        Cut <- t(Cut)
        dimnames(Cut)=list(NULL, c("CutOff", "se", "sp"))
    }
    return(Cut)
}

