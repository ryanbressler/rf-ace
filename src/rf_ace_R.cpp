#include <Rcpp.h>

RcppExport SEXP rf_ace(SEXP ns, SEXP xs) {

  int n = Rcpp::as<int>(ns);

  double x = Rcpp::as<double>(xs);

  for (int i=0; i<n; i++) {
    x=1/(1+x);
  }

  return Rcpp::wrap(x);

}
