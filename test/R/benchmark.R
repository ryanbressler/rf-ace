library(rfacer)
library(randomForest)

source("test/R/utils.R")

pMissing <- c(0.0,0.1,0.2)

errorlist <- list()

# Benchmark missing values
for ( i in 1:length(pMissing) ) {
  tmp <- benchmarkMissingValues(pMissing[[i]])
  errorlist[[i]] <- tmp$errors
} 

colors <- c("black","darkgrey","lightgrey")

pdf("errors.pdf",width=6,height=6)
par(mfcol=c(3,3))
for ( i in 1:3 ) {
tit <- paste(c("NUM, %NA=",as.character(pMissing[i]*100)),collapse='')
barplot(errorlist[[i]]$num,legend.text=FALSE,axes=TRUE,ylab="RMSE",col=colors,cex.names=0.9)
title(tit)

tit <- paste(c("NUM+TXT, %NA=",as.character(pMissing[i]*100)),collapse='')
barplot(errorlist[[i]]$txt,legend.text=FALSE,axes=TRUE,ylab="RMSE",col=colors,cex.names=0.9)
title(tit)

tit <- paste(c("NUM+CAT, %NA=",as.character(pMissing[i]*100)),collapse='')
barplot(errorlist[[i]]$cat,legend.text=FALSE,axes=TRUE,ylab="RMSE",col=colors,cex.names=0.9)
title(tit)
}
dev.off()


# Benchmark categorical splitter speed
offsets <- as.vector(c(0,2,3,4,5,7,9,15,20))
nCategories <- 3 + 3*offsets 
speeds <- data.frame(rf=vector(length=length(offsets)),rface=vector(length=length(offsets)))
RMSE <- data.frame(rf=vector(length=length(offsets)),rface=vector(length=length(offsets)))
for ( i in 1:length(offsets) ) {
  out <- benchmarkCatSplitterSpeed(offsets[i])
  speeds$rf[i] <- out$rfSpeed
  speeds$rface[i] <- out$rfaceSpeed
  RMSE$rf[i] <- out$rfRMSE
  RMSE$rface[i] <- out$rfaceRMSE
}
rownames(speeds) <- c(nCategories)
speeds <- t(speeds)
rownames(RMSE) <- c(nCategories)
RMSE <- t(RMSE)
pdf("catsplitter_speeds.pdf",width=4,height=6)
par(mfcol=c(2,1))
barplot(speeds,beside=TRUE,legend=TRUE,cex.names=0.8)
title("Categorical splitter execution time")
barplot(RMSE,beside=TRUE,legend=TRUE,cex.names=0.8)
title("RMSE as function of cardinality")
dev.off()

# Benchmark tree size
treesizes <- as.data.frame(t(read.table("tmp/treesizes.tsv")))
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

nSamples <- 1000
std <- 0.3
offset <- 0
pMissing <- 0.0

trainData <- makeData(nSamples,std,offset,pMissing)
testData <- makeData(nSamples,std,offset,pMissing)
qrface <- rface.train(trainData[c(1,2,3,5,6,7,8,9)],"N:output",nTrees=100,mTry=3,nodeSize=3,forestType="RF",quantiles=as.vector(c(0.1,0.3,0.5,0.7,0.9)))
qrfaceOut <- rface.predict(qrface,testData)

q1 <- getQuantileVector(qrfaceOut$predictions,1)
q2 <- getQuantileVector(qrfaceOut$predictions,2)
q3 <- getQuantileVector(qrfaceOut$predictions,3)
q4 <- getQuantileVector(qrfaceOut$predictions,4)
q5 <- getQuantileVector(qrfaceOut$predictions,5)

pdf("quantilepredictions.pdf")
plot(qrfaceOut$trueData,q1,col="lightgray",pch=".",cex=2)
grid()
points(qrfaceOut$trueData,q2,col="darkgray",pch=".",cex=2)
points(qrfaceOut$trueData,q3,col="black",pch=".",cex=2)
points(qrfaceOut$trueData,q4,col="darkgray",pch=".",cex=2)
points(qrfaceOut$trueData,q5,col="lightgray",pch=".",cex=2)
dev.off()

