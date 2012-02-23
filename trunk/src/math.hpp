#ifndef MATH_HPP
#define MATH_HPP

#include <vector>
#include "datadefs.hpp"
#include "errno.hpp"

using namespace std;
using datadefs::num_t;


namespace math {

  // Returns the p'th percentile of the data vector x
  num_t percentile(vector<num_t> x, const num_t p);

  /**
     Updates the mean and squared error by ADDING x_n
     to the set
     NOTE: NANs will corrupt the data
  */
  inline void incrementSquaredError(const num_t& x_n,
				    const size_t& n,
				    num_t& mu,
				    num_t& se) {
    
    // Increment sample count by one
    //++n; 
    
    // If n overflows, an error will be raised
    //if (n == 0) { 
    //  throw ERRNO_NUMERIC_OVERFLOW; 
    //}

    // Save the current mean
    num_t mu_old = mu;
    
    // Update mean
    mu += (x_n - mu_old) / n;
    
    //If there are already at least two data points, squared error can be calculated, otherwise assign se := 0.0
    if(n > 1) {
      se += (x_n - mu) * (x_n - mu_old);
    } else { 
      se = 0.0; /** Implementation note: this may spuriously invoke on
		    overflow. size_t is assumed to always be unsigned in
		    accordance with the 1999 ISO C standard (C99). */
    }
  } 

  /**
     Updates the mean and squared error by REMOVING x_n
     from the set
     NOTE: NANs will corrupt the data
  */
  /*
    inline void decrementSquaredError(const num_t& x_n,
    size_t& n,
    num_t& mu,
    num_t& se) {
    
    // Decrement sample count by one
    --n;
    
    // If n underflows, an error will be raised
    if ( n == static_cast<size_t>(-1) ) {
    throw ERRNO_NUMERIC_UNDERFLOW;
    }
    
    if(n > 1) {
    num_t mu_old = mu;
    mu -= (x_n - mu) / n;
    se -= (x_n - mu) * (x_n - mu_old);
    } else if(n == 1) {
    mu -= (x_n - mu) / n;
    se = 0.0;
    } else {
    mu = 0.0;
    se = 0.0;
    }
    
    }
  */


}
  

#endif
