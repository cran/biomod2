.testnull <-
function(object, Prev = 0.5 , dat){
  
  if( is.finite(object$deviance) & is.finite(object$null.deviance)){
    if(object$deviance != object$null.deviance){
      pred <- predict(object, dat, type="response")
    }
  }
  
  if(!exists('pred')){
    if(Prev < 0.5) pred <- rep(0, nrow(dat))
    if(Prev >= 0.5) pred <- rep(1, nrow(dat))
  }

    return(pred)
}

