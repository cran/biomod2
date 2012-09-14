.testnull <-
function(object, Prev = 0.5 , dat){

    if(object$deviance == object$null.deviance){
        if(Prev < 0.5) pred <- rep(0, nrow(dat))
        if(Prev >= 0.5) pred <- rep(1, nrow(dat))
    }
    else pred <- predict(object, dat, type="response")    
    return(pred)
}

