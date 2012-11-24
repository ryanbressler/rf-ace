rface.train <-
function(trainData, target, featureWeights = as.vector(rep(1,length(colnames(trainData)))), nTrees = 100, mTry = 10, nodeSize = 3, nMaxLeaves = 1000, nThreads = 1) {
  #if ( featureWeights == FALSE ) {
  #  featureWeights = as.vector(list(rep(1,length(colnames(trainData)))))
  #}
  predictorObj <- .Call("rfaceTrain", trainData, as.character(target), featureWeights, nTrees, mTry, nodeSize, nMaxLeaves, nThreads)
  return(predictorObj)
}
