#include <Rcpp.h>

#include <cstdlib>
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include <stdio.h>
#include <iomanip>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <set>
//#include <tuple>

#include "rf_ace.hpp"
#include "stochasticforest.hpp"
#include "treedata.hpp"
#include "datadefs.hpp"
#include "utils.hpp"
#include "options.hpp"
#include "statistics.hpp"
#include "progress.hpp"
#include "math.hpp"
#include "distributions.hpp"
#include "timer.hpp"

using namespace std;
using datadefs::num_t;

RcppExport SEXP rf_ace(SEXP ns, SEXP xs) {

  int n = Rcpp::as<int>(ns);

  double x = Rcpp::as<double>(xs);

  for (int i=0; i<n; i++) {
    x=1/(1+x);
  }

  return Rcpp::wrap(x);

}
