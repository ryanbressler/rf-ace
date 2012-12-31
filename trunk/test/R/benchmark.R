library(rfacer)

nSamples <- 1000

bag <- list(list("my","name","is","timo"),list("this","is","class2"),list("words","in","this","class3"))

classes <- sample(1:3,nSamples,replace=T)

nWordsPerSample <- sample(2:4,nSamples,replace=T)

text <- vector()

x <- 2*pi*(1:nSamples)/nSamples + rnorm(nSamples,0,0.1)
y <- sin(x) + rnorm(nSamples,0,0.1)

for ( i in 1:nSamples ) {
  c <- classes[i]
  text[i] <- paste(sample(bag[[c]],nWordsPerSample[i],replace=T),collapse=', ') 
  y[i] <- y[i] + 2.2*c  
}

n1 <- rnorm(nSamples)
n2 <- rnorm(nSamples)
# t  <- as.vector(rep(as.character("foo"),nSamples))

# Populating the data frame with the training data
data <- data.frame(y,x,n1,n2,text,as.character(classes),stringsAsFactors=FALSE)
colnames(data) <- c("N:output","N:input","N:noise1","N:noise2","T:random","C:class")

# Populating sample names
rownames(data) <- paste(c(rep("s",nSamples)),(1:nSamples),sep='')

fWeightsA <- as.vector(c(1,1,1,1,0,0))
fWeightsB <- as.vector(c(1,1,1,1,2,0))
fWeightsC <- as.vector(c(1,1,1,1,0,1))

rfmA <- rface.train(data,"N:output",featureWeights=fWeightsA,mTry=2,forestType="RF")
rfmB <- rface.train(data,"N:output",featureWeights=fWeightsB,mTry=5,forestType="RF")
rfmC <- rface.train(data,"N:output",featureWeights=fWeightsC,mTry=2,forestType="RF")

outA <- rface.predict(rfmA,data);
outB <- rface.predict(rfmB,data);
outC <- rface.predict(rfmC,data);

colors <- as.factor(data$"C:class")

dev.new()
pairs(data[1:4],col=colors)

dev.new()
par(mfcol=c(1,3))
plot(outA$predData,outA$trueData,col=colors,title="RF without textual data")
grid()
plot(outB$predData,outB$trueData,col=colors,title="RF with textual data")
grid()
plot(outC$predData,outC$trueData,col=colors,title="RF with true classes")
grid()