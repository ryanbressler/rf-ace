rface.predict <-
function(predictorObj,testData,quantiles=vector(length=0),nSamplesForQuantiles=10) {
  predictions <- .Call("rfacePredict",predictorObj,testData,quantiles,nSamplesForQuantiles);
  return(predictions)
}
