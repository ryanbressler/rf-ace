library(Rcpp)
library(rfacer)

trainData <- read.afm("test_103by300_mixed_nan_matrix.afm")

predictorObj <- rface.train(trainData,"C:class",mTry = 30, nTrees = 1000)
predictions <- rface.predict(predictorObj,trainData)
