#include <algorithm>
#include "math.hpp"

/**
   Returns the p'th percentile of the data vector x 
   NOTE: x needs to be copied due to the in-place sort operation
*/
num_t math::percentile(vector<num_t> x,
		       const num_t p) {
  
  // Initialize the percentile to NAN
  num_t prc = datadefs::NUM_NAN;
  
  // If the data vector has length 0, return
  if ( x.size() == 0 ) {
    return( prc );
  }
  
  // Sort data to increasing order
  sort(x.begin(),x.end());
  
  // Exact index without rounding
  num_t k = ( x.size() - 1 ) * p;

  // Lower bound of the index
  num_t f = floor(k);

  // Upper bound of the index
  num_t c = ceil(k);
  
  // If the upper and lower bounds are equal, 
  // we can calculate the percentile directly
  // by the index k
  if(fabs(f - c) < datadefs::EPS) {
    prc = x[static_cast<size_t>(k)];
  } else {

    // Otherwise we will interpolate linearly based on the 
    // distances from the intermediate point (k) to both 
    // bounds: ceil->k and k->floor
    num_t d0 = x[static_cast<size_t>(f)] * (c - k);
    num_t d1 = x[static_cast<size_t>(c)] * (k - f);

    // This operation equals to the weighted average,
    // which in other words is the interpolated percentile
    // we were after
    prc = d0 + d1;
  }
  
  // Finally return the calculated percentile
  return( prc );
  
}
