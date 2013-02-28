#include <algorithm>
#include <cmath>
//#include "gamma.hpp"
#include "math.hpp"

void math::transformLogistic(size_t nCategories,
			     vector<num_t>& prediction, 
			     vector<num_t>& probability) {

  //size_t nCategories = trainData_->nCategories();

  // Multiclass logistic transform of class probabilities from current probability estimates.
  assert(nCategories == prediction.size());
  vector<num_t>& expPrediction = probability; // just using the space by a different name

  // find maximum prediction
  vector<num_t>::iterator maxPrediction = max_element(prediction.begin(),prediction.end());
  // scale by maximum to prevent numerical errors

  num_t expSum = 0.0;
  size_t k;
  for (k = 0; k < nCategories; ++k) {
    expPrediction[k] = exp(prediction[k] - *maxPrediction); // scale by maximum
    expSum += expPrediction[k];
  }
  for (k = 0; k < nCategories; ++k) {
    probability[k] = expPrediction[k] / expSum;
  }
}


void math::adjustPValues(vector<num_t>& pValues, const size_t nTests) {

  num_t previousPValue = 0.0;

  for ( size_t i = 0; i < pValues.size(); ++i ) {
    
    pValues[i] *= nTests / ( i + 1 );
    
    if ( pValues[i] > 1.0 ) {
      pValues[i] = 1.0;
    }
    
    if ( pValues[i] < previousPValue ) {
      pValues[i] = previousPValue;
    } else {
      previousPValue = pValues[i];
    }
    
  }
  
}


/**
   Two-sample t-test
*/
num_t math::ttest(const vector<num_t>& x,
		  const vector<num_t>& y,
		  const bool WS) {

  // Sample mean and variance of x
  num_t mean_x = math::mean(x);
  num_t var_x = math::var(x,mean_x);
  size_t n_x = x.size();

  // If sample size is too small, we exit
  if ( n_x < 2 ) {
    return( datadefs::NUM_NAN );
  }

  // Sample mean and variance of y
  num_t mean_y = math::mean(y);
  num_t var_y = math::var(y,mean_y);
  size_t n_y = y.size();

  // If sample size is too small, we exit
  if ( n_y < 2 ) {
    return( datadefs::NUM_NAN );
  }

  // Degrees of freedom
  num_t v;

  // Standard deviation
  num_t s;

  if ( !WS ) {
    v = static_cast<num_t>( n_x + n_y - 2 );
    num_t sp = sqrt(((n_x-1) * var_x + (n_y-1) * var_y) / v);
    s = sp * sqrt(1.0 / n_x + 1.0 / n_y);
  } else {
    num_t h1 = pow(var_x / n_x + var_y / n_y,2);
    num_t h2 = pow( var_x / n_x, 2) / (n_x - 1) + pow(var_y/n_y,2)/(n_y-1);
    v = h1 / h2 ;
    s = sqrt( var_x / n_x + var_y / n_y );
  }

  // If pooled standard deviation is zero...
  if ( fabs(s) < datadefs::EPS ) {
    if ( mean_x > mean_y ) {
      return( datadefs::EPS ); // ... and x larger than y => p = EPS 
    } else if ( fabs( mean_x - mean_y ) < datadefs::EPS ) {
      return( 0.5 ); // ... and x and y almost equal => p = 0.5
    } else {
      return( 1.0 ); // ... and x smaller than y => p = 1.0
    }
  }

  // T-test statistic
  num_t tvalue = (mean_x - mean_y) / s;

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
  return( exp( lgamma(a) + lgamma(b) - lgamma(a+b) ) );
  // return( exp( LogGamma(a) + LogGamma(b) - LogGamma(a+b) ) );
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

  return( sgn*sqrt(1.0 - exp(-x2*(4.0/datadefs::NUM_PI + datadefs::A*x2) / (1+datadefs::A*x2))) );

}

num_t math::pearsonCorrelation(const vector<num_t>& x,
			       const vector<num_t>& y) {

  assert( x.size() == y.size() );
  
  size_t n = x.size();

  if ( n == 0 ) {
    return( datadefs::NUM_NAN );
  }

  num_t corr = 0.0;
  
  num_t mu_x = math::mean(x);
  num_t se_x = math::var(x,mu_x) * (n - 1);
  num_t mu_y = math::mean(y);
  num_t se_y = math::var(y,mu_y) * (n - 1);
  
  for(size_t i = 0; i < n; ++i) {
    corr += ( x[i] - mu_x ) * ( y[i] - mu_y );
  }

  return( corr / sqrt(se_x*se_y) );

}

num_t math::gamma(const vector<num_t>& x, const size_t nCategories) {
  
  size_t n = x.size();
  assert( n > 0 );
  assert( nCategories > 0 );
  
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
    return( (numerator*(nCategories - 1)) / (denominator*nCategories) );
  }
  
}


num_t math::numericalError(const vector<num_t>& x, const vector<num_t>& y) {

  assert( x.size() == y.size() );

  size_t n = x.size();

  if ( n == 0 ) {
    return(datadefs::NUM_NAN);
  }

  num_t ret = 0.0;

  for ( size_t i = 0; i < n; ++i ) {
    ret += pow( x[i] - y[i], 2 ) / n;
  }

  return( sqrt(ret) );

}


num_t math::var(const vector<num_t>& x) {

  num_t mu = math::mean(x);

  return( math::var(x,mu) );

}

num_t math::var(const vector<num_t>& x, const num_t& mu) {
  
  if ( x.size() < 2 ) {
    return( datadefs::NUM_NAN );
  }

  size_t n = x.size();

  num_t ret = 0.0;

  for(size_t i = 0; i < n; ++i) {
    ret += pow(x[i] - mu,2) / ( n - 1 );
  }

  return( ret );

}

