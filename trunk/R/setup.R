
# Only needed when running for for the first time
#install.packages("Rcpp", dependencies = TRUE);

# Load library
library("Rcpp");

# Load the dynamic library containing the C++ program
dyn.load("lib/rf_ace_R.so");

rface.train <- function(trainData, target, nTrees = 100, mTry = 10, nodeSize = 3, nMaxLeaves = 1000) {
  .Call("rfaceTrain", trainData, as.character(target), nTrees, mTry, nodeSize, nMaxLeaves);
}

rface.predict <- function(predictor,testData) {
  .Call("rfacePredict",predictor,testData);
}


