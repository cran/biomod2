Find.Optim.Stat <- function(Stat='TSS',Fit,Obs,Precision = 5, Fixed.thresh = NULL){
  if(length(unique(Obs)) == 1 | length(unique(Fit)) == 1){
    warnings("\nObserved or fited data contains only a value.. Evaluation Methods switched off\n")
    best.stat <- cutoff <- true.pos <- sensibility <- true.neg <- specificity <- NA  
  } else {
    if(Stat != 'ROC'){
      StatOptimum <- getStatOptimValue(Stat)
      if(is.null(Fixed.thresh)){ # test a range of threshold to get the one giving the best score
        mini <- max(min(quantile(Fit,0.05, na.rm=T), na.rm=T),0)
        maxi <- min(max(quantile(Fit,0.95, na.rm=T), na.rm=T),1000)
        valToTest <- unique( round(c(seq(mini,maxi,length.out=100), mini, maxi)) )
#         valToTest <- unique( c(seq(mini,maxi,by=Precision), mini, maxi) )        
      } else{
        valToTest <- Fixed.thresh
      }

      calcStat <- sapply(lapply(valToTest, function(x){return(table(Fit>x,Obs))} ), calculate.stat, stat=Stat)
      
      # rescal on 0-1 ladder.. 1 is the best
      calcStat <- 1 - abs(StatOptimum - calcStat)
      
      best.stat <- max(calcStat, na.rm=T)
      cutoff <- mean(valToTest[which(calcStat==best.stat)]) # if several values are selected

      misc <- table(Fit >= cutoff, Obs)
      misc <- .contagency.table.check(misc)
      true.pos <- misc['TRUE','1']
      true.neg <- misc['FALSE','0']
      specificity <- (true.neg * 100)/sum(misc[,'0'])
      sensibility <- (true.pos * 100)/sum(misc[,'1'])
    } else{
      best.stat <- .somers2(Fit, Obs)["C"]
      Cut <- .CutOff.Optimised(Obs, Fit)
      cutoff <- Cut[1]
      sensibility <- Cut[2]
      specificity <- Cut[3]
    }
  }
  return(cbind(best.stat,cutoff,sensibility,specificity))
}

getStatOptimValue <- function(stat){
  if(stat == 'TSS') return(1)
  if(stat == 'KAPPA') return(1)
  if(stat == 'ACCURACY') return(1)
  if(stat == 'BIAS') return(1)
  if(stat == 'POD') return(1)
  if(stat == 'FAR') return(0)
  if(stat == 'POFD') return(0)
  if(stat == 'SR') return(1)
  if(stat == 'CSI') return(1)
  if(stat == 'ETS') return(1)
  if(stat == 'HK') return(1)
  if(stat == 'HSS') return(1)
  if(stat == 'OR') return(1000000)
  if(stat == 'ORSS') return(1)
}

calculate.stat <-
function(Misc, stat='TSS')
{
  # Contagency table checking
  Misc <- .contagency.table.check(Misc)
  
  # Defining Classification index
  hits <- Misc['TRUE','1']
  misses <- Misc['FALSE','1']
  false_alarms <- Misc['TRUE','0']
  correct_negatives <- Misc['FALSE','0']
  
  total <- sum(Misc)
  forecast_1 <- sum(Misc['TRUE',])
  forecast_0 <- sum(Misc['FALSE',])
  observed_1 <- sum(Misc[,'1'])
  observed_0 <- sum(Misc[,'0'])
  
  # Calculating choosen evaluating metric
  if(stat=='TSS'){
    return( (hits/(hits+misses)) + (correct_negatives/(false_alarms+correct_negatives)) -1 )
  }
  
  if(stat=='KAPPA'){
    Po <- (1/total) * (hits + correct_negatives)
    Pe <- ((1/total)^2) * ((forecast_1 * observed_1) + (forecast_0 * observed_0))
    return( (Po - Pe) / (1-Pe) )
  }
  
  if(stat=='ACCURACY'){
    return( (hits + correct_negatives) / total)
  }
  
  if(stat=='BIAS'){
    return( (hits + false_alarms) / (hits + misses))
  }
  
  if(stat=='POD'){
    return( hits / (hits + misses))
  }
  
  if(stat=='FAR'){
    return(false_alarms/(hits+false_alarms))
  }
  
  if(stat=='POFD'){
    return(false_alarms / (correct_negatives + false_alarms))
  }
  
  if(stat=='SR'){
    return(hits / (hits + false_alarms))
  }
  
  if(stat=='CSI'){
    return(hits/(hits+misses+false_alarms))
  }
  
  if(stat=='ETS'){
    hits_rand <- ((hits+misses)*(hits+false_alarms)) / total
    return( (hits-hits_rand) / (hits+misses+false_alarms-hits_rand))
  }
  
  if(stat=='HK'){
    return((hits/(hits+misses)) - (false_alarms/(false_alarms + correct_negatives)))
  }
  
  if(stat=='HSS'){
    expected_correct_rand <- (1/total) * ( ((hits+misses)*(hits+false_alarms)) +
      ((correct_negatives + misses)*(correct_negatives+false_alarms)) )
    return((hits+correct_negatives-expected_correct_rand) / (total - expected_correct_rand))
  }
  
  if(stat=='OR'){
    return((hits*correct_negatives)/(misses*false_alarms))
  }
  
  if(stat=='ORSS'){
    return((hits*correct_negatives - misses*false_alarms) / (hits*correct_negatives + misses*false_alarms))
  }
}

.contagency.table.check <- function(Misc){
  # Contagency table checking
  if(dim(Misc)[1]==1){
    if(row.names(Misc)[1]=="FALSE"){
      Misc <- rbind(Misc, c(0,0))
      rownames(Misc) <- c('FALSE','TRUE')
    } else{
      a <- Misc
    	Misc <- c(0,0)
  		Misc <- rbind(Misc, a)
      rownames(Misc) <- c('FALSE','TRUE')
  	}
  }
  
  if(ncol(Misc) != 2 | nrow(Misc) !=2 ){
    Misc = matrix(0, ncol=2, nrow=2, dimnames=list(c('FALSE','TRUE'), c('0','1')))
  }
  
  if((sum(colnames(Misc) %in% c('FALSE','TRUE','0','1')) < 2) | (sum(rownames(Misc) %in% c('FALSE','TRUE','0','1')) < 2) ){
    stop("Unavailable contagency table given")
  }
  
  if('0' %in% rownames(Misc)) rownames(Misc)[which(rownames(Misc)=='0')] <- 'FALSE'
  if('1' %in% rownames(Misc)) rownames(Misc)[which(rownames(Misc)=='1')] <- 'TRUE'  
    
  return(Misc)
}