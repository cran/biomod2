`makeFormula` <-
function(respName, explVar, type = 'simple', interaction.level = 0)
{
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  # This function return a string in a well formated way. May be give as formula argument to a "basic"
  # statistical model.
  # Several types of models are available
  #
  # D.GEORGES 12/2011
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  
  # 0. Supported Types =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  availableTypes = c("simple", "quadratic", "polynomial", "s_smoother")
  
  # 1. Check Given Args =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if(!is.character(respName) || length(respName)!=1){
    stop("You must give a unique response variable name")
  }
  
  if(!is.data.frame(explVar) &&  !is.matrix(explVar)){
    stop("You must give explanatory variable table")
  }
  
  if(!(type %in% availableTypes)){
    stop(paste("Formuula type must be one of : ", toString(availableTypes), sep=""))
  }
  
  explVarNames <- colnames(explVar)
  if(respName %in% explVarNames){ # remove the response variable data if it's given
    explVar <- explVar[, - which(explVarNames == respName)]
    explVarNames <- colnames(explVar)
  }
  
  interaction.level <- min(interaction.level, ncol(explVar))
  
  # 2. Create the formula =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  junk <- c()  
  
  switch(EXPR=type,
         "simple" = {junk <- explVarNames},
  
         "quadratic" = { 
           for (v in 1:ncol(explVar) ){
             if(is.numeric(explVar[,v])){
#                junk <- c(junk, paste(explVarNames[v], "+I(", explVarNames[v], 
#                                       "^2)+I(",explVarNames[v],"^3)", sep="") )
                junk <- c(junk, paste(explVarNames[v], "+I(", explVarNames[v], 
                                      "^2)", sep="") )
             } else { junk <- c(junk, explVarNames[v]) }
           } },

         "polynomial" = {
           for (v in 1:ncol(explVar) ){
             if(is.numeric(explVar[,v])){
#                junk <- c(junk, paste(explVarNames[v],
#                                       "+I(", explVarNames[v],
#                                       "^2)+I(",explVarNames[v],
#                                       "^3)+poly(",explVarNames[v],
#                                       ",2)+poly(",explVarNames[v],
#                                       ",3)",sep="") )
                  junk <- c(junk, paste(explVarNames[v],
                                      "+poly(",explVarNames[v],
                                      ",3)",sep="") )               
             } else { junk <- c(junk, explVarNames[v]) }
           } },
         
         "s_smoother" = {
           for (v in 1:ncol(explVar) ){
             if(is.numeric(explVar[,v])){
                  junk <- c(junk, paste(explVarNames[v],
                                      "+s(",explVarNames[v],
                                      ")",sep="") )               
             } else { junk <- c(junk, explVarNames[v]) }
           } })
  
  # interactions
  junk.inter <- NULL
  if(interaction.level > 0){
    junk.inter <- unlist(strsplit(junk,"+",fixed=TRUE))
    eval(parse(text=paste("junk.inter <- levels(interaction(",
                          toString(rep('junk.inter', interaction.level+1)),",sep=':'))",sep="")))
    
  }
  
  junk <- gsub(', ',' + ', toString(c(junk, junk.inter)))

  # 2. Return the formula =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  return(as.formula(paste(respName," ~ ", junk, sep="")))
      
}