rface.predict <-
function(predictor,testData) {
  .Call("rfacePredict",predictor,testData);
}
