

rmse <- function(out) {
  return(sqrt(mean((out$predData-out$trueData)^2,na.rm=TRUE)))
}



sampleFakeClasses <- function(classes,offset) {

  fakeClasses <- classes

  if ( offset > 0 ) {
    for ( i in 1:length(classes) ) {
      start <- (classes[i]-1)*300 + 1
      fakeClasses[i] <- sample(start:(start+offset),1,replace=T)
    }
  }

  return(fakeClasses)
}




makeData <- function(nSamples,std,bags,nWordsMin,nWordsMax,offset,pMissing) {

  classes <- sample(1:3,nSamples,replace=T)
  fakeClasses <- sampleFakeClasses(classes,offset)

  nWordsPerSample <- sample(nWordsMin:nWordsMax,nSamples,replace=TRUE)

  text <- vector()

  v  <- seq(0,4*pi,length.out=nSamples)
  x1 <- sin(v) + rnorm(nSamples,0,std)
  x2 <- v + rnorm(nSamples,0,std)
  y  <- x1 + x2 + rnorm(nSamples,0,std)

  nNoisyVars <- 4

  for ( i in 1:nSamples ) {
    c <- classes[i]
    nWords <- nWordsPerSample[i]
    # nWords <- 10
    text[i] <- paste(sample(bags[[c]],nWords,replace=F),collapse=', ') 
    y[i] <- y[i] + 4 * pi * c  
  }

  n1 <- rnorm(nSamples)
  n1[runif(nSamples) < pMissing] <- NA
  n2 <- rnorm(nSamples)
  n2[runif(nSamples) < pMissing] <- NA
  n3 <- rnorm(nSamples)
  n3[runif(nSamples) < pMissing] <- NA
  n4 <- rnorm(nSamples)
  n4[runif(nSamples) < pMissing] <- NA
  x1[runif(nSamples) < pMissing] <- NA
  x2[runif(nSamples) < pMissing] <- NA

  # Populating the data frame with the training data
  data <- data.frame(y,x1,x2,text,as.character(fakeClasses),n1,n2,n3,n4,stringsAsFactors=FALSE)
  colnames(data) <- c("N:output","N:input1","N:input2","T:random","C:class","N:noise1","N:noise2","N:noise3","N:noise4")

  # Populating sample names
  rownames(data) <- paste(c(rep("s",nSamples)),(1:nSamples),sep='')

  return(data)

}




benchmarkMissingValues <- function(pMissing) {

nSamples <- 1000
std <- 0.4
nWordsMin <- 4
nWordsMax <- 8
offset <- 0

bags <- list(
list("buckler","shield","sword","helmet","gloves","horse","medieval","castle","joust","clown","extra","words","that","mix"),
list("swan","duck","duckling","bird","fly","pond","wings","feather","beak","legs","words","that","dont","distinguish"),
list("baby","diaper","toy","poo","pee","smile","cry","toddler","infant","play","text","that","dont","distinguish"))

trainData <- makeData(nSamples,std,bags,nWordsMin,nWordsMax,offset,pMissing)
testData <- makeData(nSamples,std,bags,nWordsMin,nWordsMax,offset,pMissing)

fWeightsA <- as.vector(c(1,1,1,0,0,1,1,1,1))
fWeightsB <- as.vector(c(1,1,1,4,0,1,1,1,1))
fWeightsC <- as.vector(c(1,1,1,0,1,1,1,1,1))

rfmA <- rface.train(trainData,"N:output",featureWeights=fWeightsA,nTrees=50,mTry=3,nodeSize=3,forestType="RF",noNABranching=TRUE)
rfmB <- rface.train(trainData,"N:output",featureWeights=fWeightsB,nTrees=50,mTry=6,nodeSize=3,forestType="RF",noNABranching=TRUE)
rfmC <- rface.train(trainData,"N:output",featureWeights=fWeightsC,nTrees=50,mTry=3,nodeSize=3,forestType="RF",noNABranching=TRUE)
rfmD <- rface.train(trainData,"N:output",featureWeights=fWeightsA,nTrees=50,mTry=3,nodeSize=3,forestType="RF",noNABranching=FALSE)
rfmE <- rface.train(trainData,"N:output",featureWeights=fWeightsB,nTrees=50,mTry=6,nodeSize=3,forestType="RF",noNABranching=FALSE)
rfmF <- rface.train(trainData,"N:output",featureWeights=fWeightsC,nTrees=50,mTry=3,nodeSize=3,forestType="RF",noNABranching=FALSE)

outA <- rface.predict(rfmA,testData)
outB <- rface.predict(rfmB,testData)
outC <- rface.predict(rfmC,testData)
outD <- rface.predict(rfmD,testData)
outE <- rface.predict(rfmE,testData)
outF <- rface.predict(rfmF,testData)

trainData$"C:class" <- as.factor(trainData$"C:class")
testData$"C:class"  <- as.factor(testData$"C:class")

imputedTrainData <- na.roughfix(trainData[c(1,2,3,5,6,7,8,9)])
imputedTestData  <- na.roughfix(testData[ c(1,2,3,5,6,7,8,9)])

rfOut1 <- randomForest(imputedTrainData[c(2,3,5,6,7,8)],y=imputedTrainData[[1]],xtest=imputedTestData[c(2,3,5,6,7,8)],ytest=imputedTestData[[1]],ntree=50,mtry=3)
rfOut2 <- randomForest(imputedTrainData[2:8],y=imputedTrainData[[1]],xtest=imputedTestData[2:8],ytest=imputedTestData[[1]],ntree=50,mtry=3)

outG <- list()
outG$trueData <- outF$trueData
outG$predData <- rfOut1$test$predicted

outH <- list()
outH$trueData <- outF$trueData
outH$predData <- rfOut2$test$predicted

colors <- testData$"C:class"

# dev.new()
pdf("scattermatrix.pdf")
pairs(testData[c(1,2,3,7)],col=colors)
dev.off()

# dev.new()
pdf("predictions.pdf",width=8,height=8)
par(mfcol=c(3,2))
plot(outA$predData,outA$trueData,col=colors,pch='.')
title("RF-ACE (binary) (A)")
lines( par()$usr[1:2], par()$usr[1:2] )
grid()
plot(outB$predData,outB$trueData,col=colors,pch='.')
title("RF-ACE (binary) (text) (B)")
lines( par()$usr[1:2], par()$usr[1:2] )
grid()
plot(outC$predData,outC$trueData,col=colors,pch='.')
title("RF-ACE (binary) (classes) (C)")
lines( par()$usr[1:2], par()$usr[1:2] )
grid()
plot(outD$predData,outD$trueData,col=colors,pch='.')
title("RF-ACE (ternary) (D)")
lines( par()$usr[1:2], par()$usr[1:2] )
grid()
plot(outE$predData,outE$trueData,col=colors,pch='.')
title("RF-ACE (ternary) (text) (E)")
lines( par()$usr[1:2], par()$usr[1:2] )
grid()
plot(outF$predData,outF$trueData,col=colors,pch='.')
title("RF-ACE (ternary) (classes) (F)")
lines( par()$usr[1:2], par()$usr[1:2] )
grid()
dev.off()

pdf("predictions_ref.pdf")
plot(outG$predData,outG$trueData,col=colors,pch='.')
title("RF (ref.)")
lines( par()$usr[1:2], par()$usr[1:2] )
grid()
dev.off()

pdf("predictions_ref2.pdf")
plot(outH$predData,outH$trueData,col=colors,pch='.')
title("RF (ref.)")
lines( par()$usr[1:2], par()$usr[1:2] )
grid()
dev.off()

#outA$method <- rep("A",nSamples)
#outB$method <- rep("B",nSamples)
#outC$method <- rep("C",nSamples)
#outD$method <- rep("D",nSamples)
#outE$method <- rep("E",nSamples)
#outF$method <- rep("F",nSamples)

#results <- data.frame(residual=c(outA$trueData-outA$predData,outB$trueData-outB$predData,outC$trueData-outC$predData,outD$trueData-outD$predData,outE$trueData-outE$predData,outF$trueData-outF$predData),method=c(outA$method,outB$method,outC$method,outD$method,outE$method,outF$method))

errors <- list()
errors$num <- c(rmse(outG),rmse(outA),rmse(outD))
names(errors$num) <- c("A","B","C")
errors$txt <- c(rmse(outG),rmse(outB),rmse(outE))
names(errors$txt) <- c("A","B","C")
errors$cat <- c(rmse(outG),rmse(outC),rmse(outF))
names(errors$cat) <- c("A","B","C")
errors$title <- paste(c("n=",as.character(nSamples), ", pMissing=",as.character(pMissing*100)),collapse='')

return(list(errors=errors,data=testData,idata=imputedTestData,rf=rfOut2,outG=outG,outH=outH))

}

benchmarkCatSplitterSpeed <- function(offset) {

nSamples <- 1000
std <- 0.3
nWordsMin <- 4
nWordsMax <- 8
pMissing <- 0.0

bags <- list(
list("buckler","shield","sword","helmet","gloves","horse","medieval","castle","joust","clown","extra","words","that","mix"),
list("swan","duck","duckling","bird","fly","pond","wings","feather","beak","legs","words","that","dont","distinguish"),
list("baby","diaper","toy","poo","pee","smile","cry","toddler","infant","play","text","that","dont","distinguish"))

trainData <- makeData(nSamples,std,bags,nWordsMin,nWordsMax,offset,pMissing)
testData <- makeData(nSamples,std,bags,nWordsMin,nWordsMax,offset,pMissing)
trainData <- trainData[c(1,5)]
testData <- testData[c(1,5)]

speed <- list()

speed$rface <- 0
for ( i in 1:10 ) {
diff <- proc.time()
rface <- rface.train(trainData,"N:output",nTrees=50,mTry=1,nodeSize=3,forestType="RF",noNABranching=FALSE)
diff <- proc.time() - diff
speed$rface <- as.matrix(speed$rface + diff)[1]
}

RMSE <- list()
rfaceOut <- rface.predict(rface,testData)
RMSE$rface <- rmse(rfaceOut)

trainData$"C:class" <- as.factor(trainData$"C:class")

#trainData <- as.matrix(na.roughfix(trainData))

speed$rf <- NA
if (offset < 10) {
speed$rf <- 0
for ( i in 1:10 ) {
diff <- proc.time()
rf <- randomForest(trainData[2],y=trainData[[1]],ntree=50,mtry=1)
diff <- proc.time() - diff
speed$rf <- as.matrix(speed$rf + diff)[1]
}
}

RMSE$rf <- NA 
if (offset < 10) {
rf <- randomForest(trainData[2],y=trainData[[1]],xtest=testData[2],ytest=testData[[1]],ntree=50,mtry=1)

rfOut <- list()
rfOut$trueData <- rfaceOut$trueData
rfOut$predData <- rf$test$predicted
RMSE$rf <- rmse(rfOut)
}

return(list(rfSpeed=speed$rf,rfaceSpeed=speed$rface,data=trainData,rfRMSE=RMSE$rf,rfaceRMSE=RMSE$rface))
}



