library(rfacer)

# Loading training data from an .afm file
trainData1 <- read.afm("test_103by300_mixed_nan_matrix.afm")

# Building predictor for "N:output"
predictorObj1 <- rface.train(trainData1,"N:output")

# 
associations1 <- rface.filter(trainData1,"N:output")
predictions1 <- rface.predict(predictorObj1,trainData1)
rface.save(predictorObj1,"foo.sf")
predictorObj2 <- rface.load("foo.sf")
predictions2 <- rface.predict(predictorObj2,trainData1)

# Generating new training data
x  <- 2*pi*(1:1000)/1000 + rnorm(1000,0,0.1)
y  <- sin(x) + rnorm(1000,0,0.1)
n1 <- rnorm(1000)
n2 <- rnorm(1000)
t  <- as.vector(rep(as.character("foo"),1000))

# Populating the data frame with the training data
trainData2 <- data.frame(y,x,n1,n2,t,stringsAsFactors=FALSE)
colnames(trainData2) <- c("N:output","N:input","N:noise1","N:noise2","T:random")

# Populating sample names
rownames(trainData2) <- paste(c(rep("s",1000)),(1:1000),sep='')

# Calculating associations for "N:output" with RF-ACE using 4 threads
associations3 <- rface.filter(trainData2,"N:output",nThreads=4)

# 
predictorObj3 <- rface.train(trainData2,"N:output")
predictions3 <- rface.predict(predictorObj3,trainData2)




