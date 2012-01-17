#ifndef MATH_HPP
#define MATH_HPP

#include <vector>
#include <algorithm>
#include "datadefs.hpp"


using namespace std;
using datadefs::num_t;


namespace math {
  
  // !! Documentation: computes the percentile for the value x, given an alpha
  // !! score and "prc." I have no idea what "prc" means, here.
  num_t percentile(vector<num_t> x,
		   const num_t alpha) {
    
    num_t prc = datadefs::NUM_NAN;
    
    sort(x.begin(),x.end());
    
    num_t k( ( x.size() - 1 ) * alpha);
    num_t f = floor(k);
    num_t c = ceil(k);
    
    if(fabs(f - c) < datadefs::EPS) {
      prc = x[static_cast<size_t>(k)];
    } else {
      num_t d0 = x[static_cast<size_t>(f)] * (c - k);
      num_t d1 = x[static_cast<size_t>(c)] * (k - f);
      prc = d0 + d1;
    }
    
    
    return( prc );
    
  }
}
  

#endif
