read.afm <- function(file) 
{

trainData <- read.table(file,head=TRUE,sep="\t",row.names=1)

trainData <- as.data.frame(t(trainData))

featureNames <- names(trainData)

for( i in 1:length(featureNames) ) {
        if ( substr(featureNames[i],1,2) != "N:" ) {
                for( j in 1:length(row.names(trainData)) ) {
                        trainData[j,i] <- as.character(trainData[j,i])
                }
        }
}
return(trainData)
}
