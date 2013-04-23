rface.predict <-
function(predictorObj,testData,quantiles=vector(length=0),nSamplesForQuantiles=10,distributions=FALSE) {
  predictions <- .Call("rfacePredict",predictorObj,testData,quantiles,nSamplesForQuantiles,distributions);
  return(predictions)
}
