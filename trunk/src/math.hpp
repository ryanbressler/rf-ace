#ifndef MATH_HPP
#define MATH_HPP

#include <vector>
#include <algorithm>
#include <map>
#include "datadefs.hpp"
#include "errno.hpp"


using namespace std;
using datadefs::num_t;


namespace math {

  /**
     Returns the p'th percentile of the data vector x
  */
  num_t percentile(vector<num_t> x, const num_t p);

  /**
     Error function
     NOTE: see http://en.wikipedia.org/wiki/Error_function
  */
  num_t erf(num_t x);

  /**
     Two-sample t-test
     NOTE: see http://en.wikipedia.org/wiki/Student's_t-test
  */
  num_t ttest(const vector<num_t>& x,
              const vector<num_t>& y);

  /**
     Regularized incomplete Beta function
     NOTE: see http://en.wikipedia.org/wiki/Beta_function
  */
  num_t regularizedIncompleteBeta(const num_t x,
				  const num_t a,
				  const num_t b);


  
  num_t pearsonCorrelation(const vector<num_t>& x,
			   const vector<num_t>& y);

  inline num_t mean(const vector<num_t>& x) {
    
    if ( x.size() == 0 ) {
      return( datadefs::NUM_NAN );
    }

    num_t mu = 0.0;

    for(size_t i = 0; i < x.size(); ++i) {
	mu += x[i];
    }

    return( mu / x.size() );
    
  }

  template<typename T>
  map<T,size_t> frequency(const vector<T>& x) {
    map<T,size_t> freq;
    for(size_t i = 0; i < x.size(); ++i) {
      if( freq.find(x[i]) == freq.end() ) {
	freq[ x[i] ] = 1;
      } else {
	++freq[ x[i] ];
      }
    } 
    return( freq );
  }

  template<typename T>
  T mode(const vector<T>& x) {
    map<T,size_t> freq = frequency<T>(x);
    typename map<T,size_t>::const_iterator maxElement( freq.begin() );
    for ( typename map<T,size_t>::const_iterator it(freq.begin()); it != freq.end(); ++it ) {
      if ( it->second > maxElement->second ) {
	maxElement = it;
      }
    }
    return( maxElement->first );
  }
  
  

  num_t gamma(const vector<num_t>& x, const size_t nCategories); 
  
  num_t squaredError(const vector<num_t>& x);
  
  num_t squaredError(const vector<num_t>& x, const num_t mu); 

  /**
     Updates the mean and squared error by ADDING x_n to the set
     NOTE: NANs will corrupt the data
  */
  inline void incrementSquaredError(const num_t& x_n,
				    const size_t& n,
				    num_t& mu,
				    num_t& se) {
    
    // Save the current mean
    num_t mu_old = mu;
    
    // Update mean
    mu += (x_n - mu_old) / n;
    
    // If there are already at least two data points, 
    // squared error can be calculated, otherwise assign se := 0.0
    if(n > 1) {
      se += (x_n - mu) * (x_n - mu_old);
    } else { 
      se = 0.0; /** Implementation note: this may spuriously invoke on
		    overflow. size_t is assumed to always be unsigned in
		    accordance with the 1999 ISO C standard (C99). */
    }
  } 


  /**
     Updates the mean and squared error by REMOVING x_n from the set
     NOTE: NANs will corrupt the data
  */
  inline void decrementSquaredError(const num_t& x_n,
				    const size_t& n,
				    num_t& mu,
				    num_t& se) {
    
    
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

  /**
     Updates the squared frequency by ADDING x_n to the set
     NOTE: NANs will corrupt the data
  */
  inline void incrementSquaredFrequency(const num_t x_n,
					map<num_t,size_t>& freq,
					size_t& sqFreq) {


    // Check if the value already exists in the frequency map
    map<num_t,size_t>::const_iterator it(freq.find(x_n));
    if(it == freq.end()) {

      // If not, squared frequency becomes updated by 1
      sqFreq += 1;
      freq[x_n] = 1;

    } else {

      // Otherwise the squared frequency becomes updated by 
      // 2*freq + 1
      sqFreq += 2*freq[x_n] + 1;
      ++freq[x_n];

    }
  }
 
  /**
     Updates the squared frequency by REMOVING x_n
     from the set
     NOTE: NANs will corrupt the data
  */ 
  inline void decrementSquaredFrequency(const num_t x_n,
                                        map<num_t,size_t>& freq,
                                        size_t& sqFreq) {
    
    assert( freq.find(x_n) != freq.end() );
    assert( freq[x_n] > 0);
    
    sqFreq -= 2*freq[x_n] - 1;
    --freq[x_n];
    
    if(freq[x_n] == 0) {
      freq.erase(x_n);
    }
  }
  
  // Calculates decrease in impurity for a numerical target
  inline num_t deltaImpurity_regr(const num_t mu_tot,
				  const size_t n_tot,
				  const num_t mu_left,
				  const size_t n_left,
				  const num_t mu_right,
				  const size_t n_right) {

    return( - mu_tot   * mu_tot
	    + mu_left  * mu_left  * n_left  / n_tot
	    + mu_right * mu_right * n_right / n_tot );
    
  }
 
  inline num_t deltaImpurity_class(const size_t sf_tot,
				   const size_t n_tot,
				   const size_t sf_left,
				   const size_t n_left,
				   const size_t sf_right,
				   const size_t n_right) {

    return( - 1.0 * sf_tot   / ( n_tot * n_tot   ) 
	    + 1.0 * sf_left  / ( n_tot * n_left  )
	    + 1.0 * sf_right / ( n_tot * n_right ) );

  }

  template <typename T>
  inline void setUnion(set<T>& baseSet, const set<T>& newSet) {

    for ( typename set<T>::const_iterator it( newSet.begin() ); it != newSet.end(); ++it ) {
      baseSet.insert(*it);
    }

  }

  
}
  
#endif
