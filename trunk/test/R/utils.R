

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




makeData <- function(nSamples,std,offset,pMissing) {

  nWordsMin <- 4
  nWordsMax <- 8

  bags <- list(
  list("buckler","shield","sword","helmet","gloves","horse","medieval","castle","joust","clown","extra","words","that","mix"),
  list("swan","duck","duckling","bird","fly","pond","wings","feather","beak","legs","words","that","dont","distinguish"),
  list("baby","diaper","toy","poo","pee","smile","cry","toddler","infant","play","text","that","dont","distinguish"))

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
  n1[runif(nSamples) < pMissing & classes == 1] <- NA
  n2 <- rnorm(nSamples)
  n2[runif(nSamples) < pMissing & classes == 1] <- NA
  n3 <- rnorm(nSamples)
  n3[runif(nSamples) < pMissing & classes == 1] <- NA
  n4 <- rnorm(nSamples)
  n4[runif(nSamples) < pMissing & classes == 1] <- NA
  x1[runif(nSamples) < pMissing & classes == 1] <- NA
  x2[runif(nSamples) < pMissing & classes == 1] <- NA
  fakeClasses[runif(nSamples) < pMissing & classes == 1] <- NA

  # Populating the data frame with the training data
  data <- data.frame(y,x1,x2,text,as.character(fakeClasses),n1,n2,n3,n4,stringsAsFactors=FALSE)
  colnames(data) <- c("N:output","N:input1","N:input2","T:random","C:class","N:noise1","N:noise2","N:noise3","N:noise4")

  # Populating sample names
  rownames(data) <- paste(c(rep("s",nSamples)),(1:nSamples),sep='')

  return(data)

}

getRFACEOutput <- function(trainData,testData,forestType,noNABranching,quantiles=vector(length(0))) {

rface <- rface.train(trainData,"N:output",nTrees=50,mTry=3,nodeSize=3,forestType=forestType,noNABranching=noNABranching)
return(rface.predict(rface,testData))

}

getQuantileVector <- function(predictions,idx) {

out <- vector(length=length(predictions))

for ( i in 1:length(predictions) ) {
out[i] <- predictions[[i]][idx]
}
return(out)
}

testCalibration <- function(rfaceOut) {

  nQuantiles <- length(rfaceOut$quantiles)
  nSamples <- length(rfaceOut$trueData)
  cal <- 1*vector(length=nQuantiles)

  for ( i in 1:nSamples ) {
    cal <- cal + 1*(rfaceOut$trueData[i] < rfaceOut$predictions[[i]])/nSamples
  }

  return(cal)
}

benchmarkMissingValues <- function(pMissing) {

offset <- 0
nSamples <- 1000
std <- 0.3

trainData <- makeData(nSamples,std,offset,pMissing)
testData <- makeData(nSamples,std,offset,pMissing)

icsNum <- as.vector(c(1,2,3,6,7,8,9))
icsNumTxt <- as.vector(c(1,2,3,4,6,7,8,9))
icsNumCat <- as.vector(c(1,2,3,5,6,7,8,9))

outA <- getRFACEOutput(trainData[icsNum],testData[icsNum],"RF",TRUE)
outB <- getRFACEOutput(trainData[icsNumTxt],testData[icsNumTxt],"RF",TRUE)
outC <- getRFACEOutput(trainData[icsNumCat],testData[icsNumCat],"RF",TRUE)
outD <- getRFACEOutput(trainData[icsNum],testData[icsNum],"RF",FALSE)
outE <- getRFACEOutput(trainData[icsNumTxt],testData[icsNumTxt],"RF",FALSE)
outF <- getRFACEOutput(trainData[icsNumCat],testData[icsNumCat],"RF",FALSE)

outG <- getRFACEOutput(trainData[icsNum],testData[icsNum],"RF",TRUE,quantiles=vector(c(0.5)))
outH <- getRFACEOutput(trainData[icsNumTxt],testData[icsNumTxt],"RF",TRUE,quantiles=vector(c(0.5)))
outI <- getRFACEOutput(trainData[icsNumCat],testData[icsNumCat],"RF",TRUE,quantiles=vector(c(0.5)))
outJ <- getRFACEOutput(trainData[icsNum],testData[icsNum],"RF",FALSE,quantiles=vector(c(0.5)))
outK <- getRFACEOutput(trainData[icsNumTxt],testData[icsNumTxt],"RF",FALSE,quantiles=vector(c(0.5)))
outL <- getRFACEOutput(trainData[icsNumCat],testData[icsNumCat],"RF",FALSE,quantiles=vector(c(0.5)))

trainData$"C:class" <- as.factor(trainData$"C:class")
testData$"C:class"  <- as.factor(testData$"C:class")

imputedTrainData <- na.roughfix(trainData[c(1,2,3,5,6,7,8,9)])
imputedTestData  <- na.roughfix(testData[ c(1,2,3,5,6,7,8,9)])

rfOut1 <- randomForest(imputedTrainData[c(2,3,5,6,7,8)],y=imputedTrainData[[1]],xtest=imputedTestData[c(2,3,5,6,7,8)],ytest=imputedTestData[[1]],ntree=50,mtry=3)
rfOut2 <- randomForest(imputedTrainData[c(2,3,4,5,6,7,8)],y=imputedTrainData[[1]],xtest=imputedTestData[c(2,3,4,5,6,7,8)],ytest=imputedTestData[[1]],ntree=50,mtry=3)

outRef1 <- list()
outRef1$trueData <- outA$trueData
outRef1$predData <- rfOut1$test$predicted

outRef2 <- list()
outRef2$trueData <- outA$trueData
outRef2$predData <- rfOut2$test$predicted

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

pdf("predictions_ref1.pdf")
plot(outRef1$predData,outRef1$trueData,col=colors,pch='.')
title("RF (ref.)")
lines( par()$usr[1:2], par()$usr[1:2] )
grid()
dev.off()

#pdf("predictions_ref2.pdf")
#plot(outH$predData,outH$trueData,col=colors,pch='.')
#title("RF (ref.)")
#lines( par()$usr[1:2], par()$usr[1:2] )
#grid()
#dev.off()

errors <- list()
errors$num <- c(rmse(outRef1),rmse(outA),rmse(outD))
names(errors$num) <- c("RF\nImputed","RF-ACE\nBinary","RF-ACE\nTernary")
errors$txt <- c(rmse(outRef1),rmse(outB),rmse(outE))
names(errors$txt) <- c("RF\nImputed","RF-ACE\nBinary","RF-ACE\nTernary")
errors$cat <- c(rmse(outRef2),rmse(outC),rmse(outF))
names(errors$cat) <- c("RF\nImputed","RF-ACE\nBinary","RF-ACE\nTernary")
errors$title <- paste(c("n=",as.character(nSamples), ", pMissing=",as.character(pMissing*100)),collapse='')

return(list(errors=errors,data=testData,idata=imputedTestData,rf=rfOut1,outG=outG,outH=outH))

}

benchmarkCatSplitterSpeed <- function(offset) {

nSamples <- 1000
std <- 0.3
pMissing <- 0.0

trainData <- makeData(nSamples,std,offset,pMissing)
testData <- makeData(nSamples,std,offset,pMissing)
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

benchmarkRFSpeeds <- function(nSamples,nIters) {
  
  speeds <- data.frame(rf=0*vector(length=length(nSamples)),qrf=0*vector(length=length(nSamples)),rface=0*vector(length=length(nSamples)))
  names(speeds) <- c("RF","QRF","RF-ACE") 
  std <- 0.3
  offset <- 0
  pMissing <- 0

  for ( iter in 1:nIters ) {
  for ( i in 1:length(nSamples) ) {

    data <- makeData(nSamples[i],std,offset,pMissing)

    diff <- proc.time()
    #rf <- randomForest(data[c(2,3,5,6,7,8,9)],y=data[[1]],ntree=50,mtry=3,nodesize=3)
    diff <- proc.time() - diff
    speeds$"RF"[i] <- as.matrix(speeds$"RF"[i]+diff)[1]

    diff <- proc.time()
    #rf <- quantregForest(data[c(2,3,5,6,7,8,9)],y=data[[1]],ntree=50,mtry=3,nodesize=3)
    diff <- proc.time() - diff
    speeds$"QRF"[i] <- as.matrix(speeds$"QRF"[i]+diff)[1]

    diff <- proc.time()
    rface <- rface.train(data[c(1,2,3,5,6,7,8,9)],"N:output",nTrees=50,mTry=3,nodeSize=3)
    diff <- proc.time() - diff
    speeds$"RF-ACE"[i] <- as.matrix(speeds$"RF-ACE"[i]+diff)[1]

  }
  }

  return(speeds)

}

