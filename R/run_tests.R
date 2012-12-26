library(rfacer)

trainData <- read.afm("test_103by300_mixed_nan_matrix.afm")

predictorObj <- rface.train(trainData,"N:output")
filterOutput <- rface.filter(trainData,"N:output")
predictions <- rface.predict(predictorObj,trainData)
rface.save(predictorObj,"foo")
predictorObj2 <- rface.load("foo")
predictions2 <- rface.predict(predictorObj2,trainData)

x <- 2*pi*(1:1000)/1000 + rnorm(1000,0,0.1); y<-sin(x); n1<-rnorm(1000);n2<-rnorm(1000); data <- data.frame("N:output"=y,"N:input"=x,"N:noise1"=n1,"N:noise2"=n2);

colnames(data) <- c("N:output","N:input","N:noise1","N:noise2")
rownames(data) <- paste(c(rep("s",1000)),(1:1000),sep='')

associations <- rface.filter(data,"N:output",nThreads=4)

predictorObj3 <- rface.train(data,"N:output")
predictions <- rface.predict(predictorObj3,data)

data["T:random"] <- rep("foo",1000)