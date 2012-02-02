#include <algorithm>
#include "math.hpp"

num_t math::percentile(vector<num_t> x,
		       const num_t alpha) {
  
  num_t prc = datadefs::NUM_NAN;
  
  if ( x.size() == 0 ) {
    return( prc );
  }
  
  sort(x.begin(),x.end());
  
  num_t k = ( x.size() - 1 ) * alpha;
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
