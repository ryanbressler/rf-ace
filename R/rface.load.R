rface.load <-
function(predictorFile, nThreads = 1) {
  predictorObj <- .Call("rfaceLoad",predictorFile,nThreads)
  return(predictorObj)
}