library(rfacer)
library(randomForest)
library(quantregForest)

source("test/R/utils.R")

pMissing <- c(0.0,0.2,0.4)

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
barplot(errorlist[[i]]$num,legend.text=FALSE,axes=TRUE,ylab="RMSE",col=colors,cex.names=0.65,ylim=c(0,9))
title(tit)

tit <- paste(c("NUM+TXT, %NA=",as.character(pMissing[i]*100)),collapse='')
barplot(errorlist[[i]]$txt,legend.text=FALSE,axes=TRUE,ylab="RMSE",col=colors,cex.names=0.65,ylim=c(0,9))
title(tit)

tit <- paste(c("NUM+CAT, %NA=",as.character(pMissing[i]*100)),collapse='')
barplot(errorlist[[i]]$cat,legend.text=FALSE,axes=TRUE,ylab="RMSE",col=colors,cex.names=0.65,ylim=c(0,9))
title(tit)
}
dev.off()


# Benchmark categorical splitter speed
offsets <- as.vector(c(0,2,3,4,5,7,9,15,20,25,30))
nCategories <- 3 + 3*offsets 
speeds <- data.frame(rf=vector(length=length(offsets)),rface=vector(length=length(offsets)))
RMSE <- data.frame(RF=vector(length=length(offsets)),"RF-ACE"=vector(length=length(offsets)))
names(RMSE) <- c("RF","RF-ACE")
for ( i in 1:length(offsets) ) {
  out <- benchmarkCatSplitterSpeed(offsets[i])
  speeds$rf[i] <- out$rfSpeed
  speeds$rface[i] <- out$rfaceSpeed
  RMSE$"RF"[i] <- out$rfRMSE
  RMSE$"RF-ACE"[i] <- out$rfaceRMSE
}
rownames(speeds) <- c(nCategories)
#speeds <- t(speeds)
rownames(RMSE) <- c(nCategories)
RMSE <- t(RMSE)
pdf("catsplitter_speeds.pdf",width=8,height=4)
par(mfcol=c(1,2))
plot(nCategories,speeds$rface,type='l',col='red',xlab='Cardinality',ylab='Runtime (s)',lwd=2.5)
lines(nCategories,speeds$rf,col='blue',lwd=2.5)
legend(10,3.5,c('RF','RF-ACE'),lty=c(1,1),lwd=c(2.5,2.5),col=c("blue","red"))
grid()
#barplot(speeds,beside=TRUE,legend=TRUE,cex.names=0.8)
#title("Categorical splitter execution time")
barplot(RMSE,beside=TRUE,legend=TRUE,cex.names=0.8,col=c("blue","red"),xlab="Cardinality",ylab="RMSE")
#plot(nCategories,RMSE$rf)
#points(nCategories,RMSE$rface)
#title("RMSE as function of cardinality")
box()
grid()
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
pMissing <- 0.2

trainData <- makeData(nSamples,std,offset,pMissing)
testData <- makeData(nSamples,std,offset,pMissing)
qrface <- rface.train(trainData[c(1,2,3,5,6,7,8,9)],"N:output",nTrees=100,mTry=3,nodeSize=10)
qrfaceOut <- rface.predict(qrface,testData,quantiles=as.vector(seq(0.01,0.90,0.05)),nSamplesForQuantiles=50)

cal <- testCalibration(qrfaceOut)

pdf("quantilecalibration.pdf")
plot(qrfaceOut$quantiles,cal,pch=".",cex=4)
lines( par()$usr[1:2], par()$usr[1:2] )
grid()
dev.off()

speeds <- benchmarkRFSpeeds()

