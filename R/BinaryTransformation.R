setGeneric("BinaryTransformation",
           function(data, threshold){
             standardGeneric("BinaryTransformation")
           })

setMethod('BinaryTransformation', signature(data='data.frame'), 
  function(data, threshold)
  {
  	FUN2 <- function(x,y){
  		moa <- apply((x>y),2,as.integer)
  		if(ncol(moa)==1) return(moa[,1])
  		else return(moa)
  	}
  	return(sweep(data.matrix(data), 2, threshold, FUN2))
  })

setMethod('BinaryTransformation', signature(data='matrix'), 
  function(data, threshold)
  {
    data <- as.data.frame(data)
    return(BinaryTransformation(data, threshold))
  })

setMethod('BinaryTransformation', signature(data='numeric'), 
  function(data, threshold)
  {
    data <- as.data.frame(data)
    return(BinaryTransformation(data, threshold))
  })

setMethod('BinaryTransformation', signature(data='RasterLayer'), 
  function(data, threshold)
  {
    return(reclassify(data,c(-Inf,threshold,0, threshold,+Inf,1)))
  })

setMethod('BinaryTransformation', signature(data='RasterStack'), 
  function(data, threshold)
  {
    if(length(threshold) == 1){
      threshold <- rep(threshold, raster:::nlayers(data))
    }
    StkTmp <- raster:::stack()
    for(i in 1:raster:::nlayers(data)){
      StkTmp <- raster:::addLayer(StkTmp, BinaryTransformation(raster:::subset(data,i,drop=TRUE), threshold[i]))
    }
    return(StkTmp)
  })
          
setMethod('BinaryTransformation', signature(data='RasterBrick'), 
  function(data, threshold)
  {
    data <- raster:::stack(data)
    return(BinaryTransformation(data, threshold))
  })