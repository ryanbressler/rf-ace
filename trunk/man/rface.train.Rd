\name{rface.train}
\alias{rface.train}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 foo bar
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rface.train(trainData, target, nTrees = 100, mTry = 10, nodeSize = 3, nMaxLeaves = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trainData}{
%%     ~~Describe \code{trainData} here~~
}
  \item{target}{
%%     ~~Describe \code{target} here~~
}
  \item{nTrees}{
%%     ~~Describe \code{nTrees} here~~
}
  \item{mTry}{
%%     ~~Describe \code{mTry} here~~
}
  \item{nodeSize}{
%%     ~~Describe \code{nodeSize} here~~
}
  \item{nMaxLeaves}{
%%     ~~Describe \code{nMaxLeaves} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (trainData, target, nTrees = 100, mTry = 10, nodeSize = 3, 
    nMaxLeaves = 1000) 
{
    .Call("rfaceTrain", trainData, as.character(target), nTrees, 
        mTry, nodeSize, nMaxLeaves)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
