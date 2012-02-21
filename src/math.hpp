#ifndef MATH_HPP
#define MATH_HPP

#include <vector>
#include "datadefs.hpp"

using namespace std;
using datadefs::num_t;


namespace math {

  // Returns the p'th percentile of the data vector x
  num_t percentile(vector<num_t> x, const num_t p);
    
}
  

#endif
