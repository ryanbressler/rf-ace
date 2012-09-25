
# Helper function to check if a package is installed
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 

# Only needed when running for for the first time
if ( ! is.installed("Rcpp") ) {
  install.packages("Rcpp", dependencies = TRUE);
}

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


