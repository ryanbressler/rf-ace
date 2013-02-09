rface.train <-
function(trainData, target, featureWeights = as.vector(rep(1,length(colnames(trainData)))), forestType = "RF", nTrees = 100, mTry = 10, nodeSize = 3, nMaxLeaves = 1000, shrinkage = 0.01, noNABranching = FALSE, quantiles = vector(length=0), nThreads = 1) {
  predictorObj <- .Call("rfaceTrain", trainData, as.character(target), featureWeights, as.character(forestType), nTrees, mTry, nodeSize, nMaxLeaves, shrinkage, noNABranching, quantiles, nThreads)
  return(predictorObj)
}
