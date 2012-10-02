rface.train <-
function(trainData, target, nTrees = 100, mTry = 10, nodeSize = 3, nMaxLeaves = 1000) {
  .Call("rfaceTrain", trainData, as.character(target), nTrees, mTry, nodeSize, nMaxLeaves);
}
