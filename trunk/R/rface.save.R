rface.save <-
function(predictorObj,fileName) {
  .Call("rfaceSave",predictorObj,fileName)
}