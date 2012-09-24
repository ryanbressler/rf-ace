
# Only needed when running for for the first time
#install.packages("Rcpp", dependencies = TRUE);

# Load library
library("Rcpp");

# Load the dynamic library containing the C++ program
dyn.load("src/rf_ace.so");

rf_ace <- function(a,b) {
  .Call("rf_ace",a,b);
}


