#include <algorithm>
#include "gamma.hpp"
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

/**
   Two-sample t-test
*/
num_t math::ttest(const vector<num_t>& x,
		  const vector<num_t>& y) {

  // Sample mean and variance of x
  num_t mean_x = 0;
  num_t var_x = 0;
  size_t n_x = x.size();
  math::squaredError(x, mean_x, var_x);

  // If sample size is too small, we exit
  if ( n_x < 2 ) {
    return( datadefs::NUM_NAN );
  }

  // Normalize squared error to get unbiased variance estimate
  var_x /= (n_x - 1);

  // Sample mean and variance of y
  num_t mean_y = 0;
  num_t var_y = 0;
  size_t n_y = y.size();
  math::squaredError(y, mean_y, var_y);

  // If sample size is too small, we exit
  if ( n_y < 2 ) {
    return( datadefs::NUM_NAN );
  }

  // Normalize squared error to get unbiased variance estimate
  var_y /= (n_y - 1);

  // Degrees of freedom
  size_t v = n_x + n_y - 2;

  // Pooled standard deviation
  num_t sp = sqrt(((n_x-1) * var_x + (n_y-1) * var_y) / v);

  // If pooled standard deviation is zero...
  if ( fabs(sp) < datadefs::EPS ) {
    if ( mean_x > mean_y ) {
      return( datadefs::EPS ); // ... and x larger than y => p = EPS 
    } else if ( fabs( mean_x - mean_y ) < datadefs::EPS ) {
      return( 0.5 ); // ... and x and y almost equal => p = 0.5
    } else {
      return( 1.0 ); // ... and x smaller than y => p = 1.0
    }
  }

  // T-test statistic
  num_t tvalue = (mean_x - mean_y) / (sp * sqrt(1.0 / n_x + 1.0 / n_y));

  // Transformed t-test statistic
  num_t ttrans = v / ( pow(tvalue,2) + v ); 
  
  // This variable will store the integral of the tail of the t-distribution
  num_t integral;

  // When ttrans > 0.9, we need to recast the integration in order to retain
  // accuracy. In other words we make use of the following identity:
  //
  // I(x,a,b) = 1 - I(1-x,b,a)
  if ( ttrans > 0.9 ) {

    // Calculate I(x,a,b) as 1 - I(1-x,b,a)
    integral = 1 - math::regularizedIncompleteBeta(1 - ttrans, 0.5, v/2);

  } else {

    // Calculate I(x,a,b) directly
    integral = math::regularizedIncompleteBeta(ttrans, v/2, 0.5);
  }

  // We need to be careful about which way to calculate the integral so that it represents 
  // the tail of the t-distribution. The sign of the tvalue hints which way to integrate
  if ( tvalue > 0.0 ) {
    return( integral / 2 );
  } else {
    return( 1 - integral / 2 );
  }

}

/**
   Odd factors for the infinite continued fraction representation of the 
   regularized incomplete beta function
*/
num_t dO(const num_t m, 
	 const num_t x, 
	 const num_t a, 
	 const num_t b) {
  return( -1.0*(a+m)*(a+b+m)*x / ( (a+2*m)*(a+2*m+1) ) );
}

/**
   Even factors for the infinite continued fraction representation of the
   regularized incomplete beta function
*/
num_t dE(const num_t m, 
	 const num_t x, 
	 const num_t a, 
	 const num_t b) {
  return( m*(b-m)*x / ((a+2*m-1)*(a+2*m)) );
}

/**
   Beta function, implemented as function of log-gamma functions implemented 
   in "gamma.hpp"
*/
num_t beta(const num_t a, const num_t b) {
  return( exp( LogGamma(a) + LogGamma(b) - LogGamma(a+b) ) );
}

// http://en.wikipedia.org/wiki/Beta_function
// http://en.wikipedia.org/wiki/Student's_t-distribution
// http://www.boost.org/doc/libs/1_38_0/libs/math/doc/sf_and_dist/html/math_toolkit/special/sf_beta/ibeta_function.html
// http://www.mpi-hd.mpg.de/astrophysik/HEA/internal/Numerical_Recipes/f6-4.pdf
num_t math::regularizedIncompleteBeta(const num_t x, 
				      const num_t a,
				      const num_t b) {

  // Number of factors in the infinite continued fraction representation
  size_t i = 50;

  num_t continuedFraction = 1; 

  // Accumulate the continued fraction
  while ( i >= 1 ) {
    num_t m = static_cast<num_t>(i);
    continuedFraction = 1 + dE(m,x,a,b) / ( 1 + dO(m,x,a,b) / continuedFraction );
    --i;
  }
  
  return( pow(x,a)*pow(1-x,b) / ( a * beta(a,b) * ( 1 + dO(0,x,a,b) / continuedFraction ) ) );

}

num_t math::erf(const num_t x) {

  num_t x2 = x*x;

  num_t sgn;
  if(x < 0.0) {
    sgn = -1.0;
  } else {
    sgn = 1.0;
  }

  return( sgn*sqrt(1.0 - exp(-x2*(4.0/datadefs::PI+datadefs::A*x2) / (1+datadefs::A*x2))) );

}

num_t math::pearsonCorrelation(const vector<num_t>& x,
			       const vector<num_t>& y) {

  assert( x.size() == y.size() );
  
  size_t n = x.size();

  if ( n == 0 ) {
    return( datadefs::NUM_NAN );
  }

  num_t corr = 0.0;
  
  num_t mu_x,se_x,mu_y,se_y;
  
  math::squaredError(x,mu_x,se_x);
  math::squaredError(y,mu_y,se_y);

  for(size_t i = 0; i < n; ++i) {
    corr += ( x[i] - mu_x ) * ( y[i] - mu_y );
  }

  return( corr / sqrt(se_x*se_y) );

}

num_t math::gamma(const vector<num_t>& x, const size_t nCategories) {
  
  size_t n = x.size();
  assert(n > 0);
  
  num_t numerator = 0.0;
  num_t denominator = 0.0;
  
  for (size_t i = 0; i < n; ++i) {
    num_t abs_data_i = fabs( x[i] );
    denominator += abs_data_i * (1.0 - abs_data_i);
    numerator   += x[i];
  }
  
  if ( fabs(denominator) <= datadefs::EPS ) {
    return( datadefs::LOG_OF_MAX_NUM * numerator );
  } else {
    return( 1.0*(nCategories - 1)*numerator / (1.0*nCategories*denominator) );
  }
  
}

