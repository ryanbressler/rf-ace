\name{rface.filter}
\alias{rface.filter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Apply feature selection with RF-ACE.
}
\description{
 Apply feature selection with RF-ACE.
}
\usage{
 associations <- rface.filter(filterData, target, featureWeights = as.vector(rep(1,length(colnames(filterData)))), nTrees = 100, mTry = 10, nodeSize = 3, nMaxLeaves = 1000, nThreads = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filterData}{
  A data.frame storing the data for feature selection. filterData[featureName] stores a vector of strings (categorical feature) or floats (numerical feature).
}
  \item{target}{
  An integer or string, pointing to a feature in filterData, i.e. filterData[target].
}
  \item{featureWeights}{
  A vector of nonnegative weights for the features; affects the sampling distribution. By default all features get weight 1 corresponding uniform sampling.
}
  \item{nTrees}{  
  Number of trees in the forest. Default 100.
}
  \item{mTry}{
  Number of randomly sampled candidate features per split. Default 10. 
}
  \item{nodeSize}{
  Minimum number of train samples per node. Default 3.
}
  \item{nMaxLeaves}{
  Maximum number of leaves per tree. Default 1000.
}
  \item{nThreads}{
  Number of CPU threads to train the model with. Default 1.
}
}
\details{

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  http://code.google.com/p/rf-ace
}
\author{
Timo Erkkila
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{ \link{read.afm}, \link{rface.train}, \link{rface.predict}, \link{rface.save}, \link{rface.load} }
}
\examples{

  afmFile <- "test_103by300_mixed_nan_matrix.afm";
  target <- "N:output";
  
  nTrees <- 100;
  mTry <- 30;
  
  nThreads <- 4;

  filterData <- read.afm(afmFile);

  associations <- rface.filter(filterData, target, nTrees = nTrees, mTry = mTry, nThreads = nThreads);

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ read.afm }
\keyword{ rface.filter }
