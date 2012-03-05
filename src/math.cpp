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

  if ( n_x < 2 ) {
    return( datadefs::NUM_NAN );
  }

  var_x /= (n_x - 1);

  // Sample mean and variance of y
  num_t mean_y = 0;
  num_t var_y = 0;
  size_t n_y = y.size();
  math::squaredError(y, mean_y, var_y);

  if ( n_y < 2 ) {
    return( datadefs::NUM_NAN );
  }

  var_y /= (n_y - 1);

  size_t v;
  num_t sp,tvalue,ttrans;

  v = n_x + n_y - 2;
  sp = sqrt(((n_x-1) * var_x + (n_y-1) * var_y) / v);
  tvalue = (mean_x - mean_y) / (sp * sqrt(1.0 / n_x + 1.0 / n_y));

  if ( tvalue > 100 ) {
    return( 0.0 );
  }

  if ( tvalue < -100 ) {
    return( 1.0 );
  }

  if ( fabs(tvalue) < datadefs::EPS ) {
    return( 0.5 );
  }

  ttrans = v / ( pow(tvalue,2) + v ); 
  
  //ttrans = (tvalue+sqrt(pow(tvalue,2) + v)) / (2 * sqrt(pow(tvalue,2) + v));

  // If tvalue is negligible to v, ttrans will become close to one, which would
  // yield inaccurate estimates for the regularized incomplete beta function
  if ( ttrans < 0.05 ) {
    cout << "woops" << endl;
    return( 0.5 );
  }

  // Approximate regularized incomplete beta function
  num_t integral = math::regularizedIncompleteBeta(ttrans, v/2, 0.5);
  
  if ( tvalue > 0.0 ) {
    return( integral / 2 );
  } else {
    return( 1 - integral / 2 );
  }

}

num_t dO(const num_t m, 
	 const num_t x, 
	 const num_t a, 
	 const num_t b) {
  return( -1.0*(a+m)*(a+b+m)*x / ( (a+2*m)*(a+2*m+1) ) );
}

num_t dE(const num_t m, 
	 const num_t x, 
	 const num_t a, 
	 const num_t b) {
  return( m*(b-m)*x / ((a+2*m-1)*(a+2*m)) );
}

num_t beta(const num_t a, const num_t b) {
  return( exp( lgamma(a) + lgamma(b) - lgamma(a+b) ) );
}

// http://en.wikipedia.org/wiki/Beta_function
// http://en.wikipedia.org/wiki/Student's_t-distribution
// http://www.boost.org/doc/libs/1_38_0/libs/math/doc/sf_and_dist/html/math_toolkit/special/sf_beta/ibeta_function.html
// http://www.mpi-hd.mpg.de/astrophysik/HEA/internal/Numerical_Recipes/f6-4.pdf
num_t math::regularizedIncompleteBeta(const num_t x, 
				      const num_t a,
				      const num_t b) {

  size_t i = 100;

  num_t foo = 1; // + dE(m,x,a,b) / ( 1 + dO(m,x,a,b) );

  while ( i >= 1 ) {
    num_t m = static_cast<num_t>(i);
    foo = 1 + dE(m,x,a,b) / ( 1 + dO(m,x,a,b) / foo );
    --i;
  }
  
  return( pow(x,a)*pow(1-x,b)/a/beta(a,b)/(1+dO(0,x,a,b)/foo) );

  //return( pow(x,a)*pow(1-x,b)/a/beta(a,b)*1/(1+ dO(0,x,a,b)/(1+ dE(1,x,a,b)/(1+ dO(1,x,a,b)/(1+ dE(2,x,a,b)/(1+ dO(2,x,a,b)))))) );

}

/*
  num_t math::regularizedIncompleteBeta(const num_t x,
  const size_t a,
  const size_t b) {
  
  num_t ibval = 0.0;
  
  size_t ab = a + b;
  
  num_t jfac = 0.0;
  for(size_t i = 2; i < a; ++i) {
  jfac += log(static_cast<num_t>(i));
  }
  
  num_t kfac = 0.0;
  for(size_t i = a+1; i < ab; ++i) {
  kfac += log(static_cast<num_t>(i));
  }
  
  for(size_t i = a; i < ab; ++i) {
  jfac += log(static_cast<num_t>(i));
  kfac += log(static_cast<num_t>(ab - i));
  num_t temp = kfac - jfac + i*log(x) + (ab-1-i)*log(1-x);
  ibval += exp(temp);
  
  //cout << jfac << "\t" << kfac << "\t" << ibval << endl;
  }
  
  if(ibval > 1.0) {
  ibval = 1.0;
  }
  
  return( ibval );
  
  }
*/

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

