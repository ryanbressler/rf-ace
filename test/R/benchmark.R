library(rfacer)

nSamples <- 10000

bag <- list(
list("my","name","is","timo","and","this","is","bag","of","words"),
list("in","this","bag","there","are","some","random","words"),
list("but","the","bag","of","words","can","be","much","larger","like","this","one","with","some","extra","random","words"))

classes <- sample(1:3,nSamples,replace=T)

nWordsPerSample <- sample(2:8,nSamples,replace=T)

text <- vector()

x <- 2*pi*(1:nSamples)/nSamples + rnorm(nSamples,0,0.1)
y <- sin(x) + rnorm(nSamples,0,0.1)

for ( i in 1:nSamples ) {
  c <- classes[i]
  text[i] <- paste(sample(bag[[c]],nWordsPerSample[i],replace=T),collapse=', ') 
  y[i] <- y[i] + 2.2 * c  
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

rfmA <- rface.train(data,"N:output",featureWeights=fWeightsA,mTry=2,nodeSize=10,forestType="RF")
rfmB <- rface.train(data,"N:output",featureWeights=fWeightsB,mTry=4,nodeSize=10,forestType="RF")
rfmC <- rface.train(data,"N:output",featureWeights=fWeightsC,mTry=2,nodeSize=10,forestType="RF")

outA <- rface.predict(rfmA,data);
outB <- rface.predict(rfmB,data);
outC <- rface.predict(rfmC,data);

colors <- as.factor(data$"C:class")

data$"C:class" <- as.factor(data$"C:class")

pdf("scattermatrix.pdf")
pairs(data[c(1,2,3,4,6)],col=colors)
dev.off()

pdf("predictions.pdf")
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
dev.off()

