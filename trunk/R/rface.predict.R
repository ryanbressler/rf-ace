rface.predict <-
function(predictorObj,testData) {
  predictions <- .Call("rfacePredict",predictorObj,testData);
  return(predictions)
}
