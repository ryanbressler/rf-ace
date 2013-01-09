library(rfacer)
library(randomForest)

generateBenchmark <- function(nSamples,std,bags,nWordsMin,nWordsMax,pMissing) {

classes <- sample(1:3,nSamples,replace=T)

nWordsPerSample <- sample(nWordsMin:nWordsMax,nSamples,replace=TRUE)

text <- vector()

v  <- seq(0,4*pi,length.out=nSamples)
x1 <- sin(v) + rnorm(nSamples,0,std)
x2 <- v + rnorm(nSamples,0,std)
y  <- x1 + x2 + rnorm(nSamples,0,std)

nNoisyVars <- 4
# pMissing <- 0.05

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
x1[runif(nSamples) < pMissing] <- NA
x2[runif(nSamples) < pMissing] <- NA

# Populating the data frame with the training data
data <- data.frame(y,x1,x2,n1,n2,text,as.character(classes),stringsAsFactors=FALSE)
colnames(data) <- c("N:output","N:input1","N:input2","N:noise1","N:noise2","T:random","C:class")

# Populating sample names
rownames(data) <- paste(c(rep("s",nSamples)),(1:nSamples),sep='')

fWeightsA <- as.vector(c(1,1,1,1,1,0,0))
fWeightsB <- as.vector(c(1,1,1,1,1,4,0))
fWeightsC <- as.vector(c(1,1,1,1,1,0,1))

rfmA <- rface.train(data,"N:output",featureWeights=fWeightsA,nTrees=50,mTry=2,nodeSize=3,forestType="RF",noNABranching=TRUE)
rfmB <- rface.train(data,"N:output",featureWeights=fWeightsB,nTrees=50,mTry=6,nodeSize=3,forestType="RF",noNABranching=TRUE)
rfmC <- rface.train(data,"N:output",featureWeights=fWeightsC,nTrees=50,mTry=2,nodeSize=3,forestType="RF",noNABranching=TRUE)
rfmD <- rface.train(data,"N:output",featureWeights=fWeightsA,nTrees=50,mTry=2,nodeSize=3,forestType="RF",noNABranching=FALSE)
rfmE <- rface.train(data,"N:output",featureWeights=fWeightsB,nTrees=50,mTry=6,nodeSize=3,forestType="RF",noNABranching=FALSE)
rfmF <- rface.train(data,"N:output",featureWeights=fWeightsC,nTrees=50,mTry=2,nodeSize=3,forestType="RF",noNABranching=FALSE)

idata <- as.matrix(na.roughfix(data[c(1,2,3,4,5)]))

outA <- rface.predict(rfmA,data)
outB <- rface.predict(rfmB,data)
outC <- rface.predict(rfmC,data)
outD <- rface.predict(rfmD,data)
outE <- rface.predict(rfmE,data)
outF <- rface.predict(rfmF,data)

rfOut <- randomForest(idata[,2:5],y=idata[,1],xtest=idata[,2:5],ytest=idata[,1],ntree=50,mtry=2)

outG <- list()
outG$trueData <- outF$trueData
outG$predData <- rfOut$predicted

colors <- as.factor(data$"C:class")

data$"C:class" <- as.factor(data$"C:class")

# dev.new()
pdf("scattermatrix.pdf")
pairs(data[c(1,2,3,7)],col=colors)
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

outA$method <- rep("A",nSamples)
outB$method <- rep("B",nSamples)
outC$method <- rep("C",nSamples)
outD$method <- rep("D",nSamples)
outE$method <- rep("E",nSamples)
outF$method <- rep("F",nSamples)

results <- data.frame(residual=c(outA$trueData-outA$predData,outB$trueData-outB$predData,outC$trueData-outC$predData,outD$trueData-outD$predData,outE$trueData-outE$predData,outF$trueData-outF$predData),method=c(outA$method,outB$method,outC$method,outD$method,outE$method,outF$method))

boxplot(residual~method,data=results)
title("Model residuals")
dev.off()

rmse <- function(out) {
  return(sqrt(mean((out$predData-out$trueData)^2,na.rm=TRUE)))
}

titleStr <- paste(c("Simulated data, ",as.character(nSamples)," samples, ",as.character(pMissing*100),"% missing values"),collapse='')

pdfStr <- paste(c("errors_",nSamples,"_",pMissing*100,".pdf"),collapse='')

pdf(pdfStr,width=9,height=4)
errors <- c(rmse(outG),rmse(outA),rmse(outD),rmse(outB),rmse(outE),rmse(outC),rmse(outF))
names(errors) <- c("RF (ref.)\n","RF-ACE\nBinary","RF-ACE\nTernary","RF-ACE\nBinary","RF-ACE\nTernary","RF-ACE\nBinary","RF-ACE\nTernary")
colors <- c("black","blue","blue","green","green","magenta","magenta")
barplot(errors,legend.text=FALSE,axes=TRUE,xlab="Models",ylab="Root Mean Squared Error",col=colors,cex.names=0.9,ylim=c(0,13))
title(titleStr)
legend("topright",c("Num.","Num.","Num. + Text","Num. + Cat."),col=c("black","blue","green","magenta"),pch=c(15,15,15,15))
dev.off()

return(list(outG=outG,data=data,idata=idata))

}

nSamples <- 1000
std <- 0.1
pMissing <- 0.05
nWordsMin <- 4
nWordsMax <- 8

bags <- list(
list("buckler","shield","sword","helmet","gloves","horse","medieval","castle","joust","clown","extra","words","that","mix"),
list("swan","duck","duckling","bird","fly","pond","wings","feather","beak","legs","words","that","dont","distinguish"),
list("baby","diaper","toy","poo","pee","smile","cry","toddler","infant","play","text","that","dont","distinguish"))

generateBenchmark(nSamples,0.3,bags,nWordsMin,nWordsMax,0.00)
generateBenchmark(nSamples,0.3,bags,nWordsMin,nWordsMax,0.10)
generateBenchmark(nSamples,0.3,bags,nWordsMin,nWordsMax,0.20)

treesizes <- as.data.frame(t(read.table("tmp/treesizes.tsv")))
# treesizes <- data.frame("0%"=tmp[1,],"10%"=tmp[2,],"20%"=tmp[3,],"30%"=tmp[4,])
names(treesizes) <- c("0%","10%","20%","30%","0%","10%","20%","30%")

pdf("treesizes.pdf",width=8,height=4)
par(mfcol=c(1,2))
boxplot(treesizes[1:4],ylim=c(15,150),xlab="% of missing values",ylab="Tree size (nodes)")
title("RF-ACE, binary splits")
grid()
boxplot(treesizes[5:8],ylim=c(15,150),xlab="% of missing values",ylab="Tree size (nodes)")
title("RF-ACE, ternary splits")
grid()
dev.off()



