rface.predict <-
function(predictorObj,testData,nThreads = 1) {
  predictions <- .Call("rfacePredict",predictorObj,testData,nThreads);
  return(predictions)
}
