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
for ( i in 1:length(offsets) ) {
  out <- benchmarkCatSplitterSpeed(offsets[i])
  speeds$rf[i] <- out$rf
  speeds$rface[i] <- out$rface
}
pdf("catsplitter_speeds.pdf",width=4,height=4)
plot(nCategories,speeds$rface,type="l",col="red")
lines(nCategories,speeds$rf,col="blue")
grid()
title("Categorical splitter execution time")
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


trainData <- read.afm("test_103by300_mixed_nan_matrix.afm")
rface <- rface.train(trainData,"N:output",nTrees=50,mTry=3,nodeSize=3,forestType="RF",quantiles=vector(0.1,0.3,0.5,0.7,0.9))


