rface.train <-
function(trainData, target, featureWeights = vector(length=0), forestType = "RF", nTrees = 100, mTry = 10, nodeSize = 3, nMaxLeaves = 1000, shrinkage = 0.01, noNABranching = FALSE, nThreads = 1) {
  predictorObj <- .Call("rfaceTrain", trainData, as.character(target), featureWeights, as.character(forestType), nTrees, mTry, nodeSize, nMaxLeaves, shrinkage, noNABranching, nThreads)
  return(predictorObj)
}
