rface.filter <-
function(filterData, target, featureWeights = vector(length=0), nTrees = 100, mTry = 10, nodeSize = 3, nMaxLeaves = 1000, nThreads = 1) {
  filterOutput <- .Call("rfaceFilter", filterData, as.character(target), featureWeights, nTrees, mTry, nodeSize, nMaxLeaves, nThreads)
  return(filterOutput)
}
