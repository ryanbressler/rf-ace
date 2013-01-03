library(rfacer)
library(randomForest)

nSamples <- 1000

bag_easy <- list(
list("buckler"),
list("swan"),
list("baby"))

bag <- list(
list("buckler","shield","sword","helmet","gloves","horse","medieval","castle","joust","clown"),
list("swan","duck","duckling","bird","fly","pond","wings","feather","beak","legs"),
list("baby","diaper","toy","poo","pee","smile","cry","toddler","infant","play"))

classes <- sample(1:3,nSamples,replace=T)

nWordsPerSample <- sample(10:10,nSamples,replace=TRUE)

text <- vector()

v <- seq(0,4*pi,length.out=nSamples)
x1 <- sin(v) + rnorm(nSamples,0,0.1)
x2 <- v + rnorm(nSamples,0,0.1)
y  <- x1 + x2 + rnorm(nSamples,0,0.1)

nVars <- 4
pMissing <- 0.05

for ( i in 1:nSamples ) {
  c <- classes[i]
  # nWords <- nWordsPerSample[i]
  nWords <- 10
  text[i] <- paste(sample(bag[[c]],nWords,replace=F),collapse=', ') 
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

rfmA <- rface.train(data,"N:output",featureWeights=fWeightsA,mTry=2,nodeSize=3,forestType="RF")
rfmB <- rface.train(data,"N:output",featureWeights=fWeightsB,mTry=6,nodeSize=3,forestType="RF")
rfmC <- rface.train(data,"N:output",featureWeights=fWeightsC,mTry=2,nodeSize=3,forestType="RF")

outA <- rface.predict(rfmA,data);
outB <- rface.predict(rfmB,data);
outC <- rface.predict(rfmC,data);

colors <- as.factor(data$"C:class")

data$"C:class" <- as.factor(data$"C:class")

dev.new()
#pdf("scattermatrix.pdf")
pairs(data[c(1,2,3,7)],col=colors)
# dev.off()

dev.new()
#pdf("predictions.pdf")
par(mfcol=c(2,2))
plot(outA$predData,outA$trueData,col=colors,pch='.')
title("RF without textual data (A)")
lines( par()$usr[1:2], par()$usr[1:2] )
grid()
plot(outB$predData,outB$trueData,col=colors,pch='.')
title("RF with textual data (B)")
lines( par()$usr[1:2], par()$usr[1:2] )
grid()
plot(outC$predData,outC$trueData,col=colors,pch='.')
title("RF with true classes (C)")
lines( par()$usr[1:2], par()$usr[1:2] )
grid()

outA$method <- rep("A",nSamples)
outB$method <- rep("B",nSamples)
outC$method <- rep("C",nSamples)

results <- data.frame(residual=c(outA$trueData-outA$predData,outB$trueData-outB$predData,outC$trueData-outC$predData),method=c(outA$method,outB$method,outC$method))

boxplot(residual~method,data=results)
title("Model residuals")
# dev.off()

