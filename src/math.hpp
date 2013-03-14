#ifndef MATH_HPP
#define MATH_HPP

#include <vector>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <utility>
#include "datadefs.hpp"
#include "errno.hpp"


using namespace std;
using datadefs::num_t;


namespace math {

  /**
     Returns the p'th percentile of the data vector x
  */
  template<class T>
  T percentile(vector<T> x, const num_t p) {

    // If the data vector has length 0, return
    if ( x.size() == 0 ) {
      cerr << "math::percentile() -- cannot compute with vector of length 0!" << endl;
      exit(1);
    }

    T prc;

    // Sort data to increasing order
    sort(x.begin(),x.end());

    // Exact index without rounding
    T k = ( x.size() - 1 ) * p;

    // Lower bound of the index
    T f = floor(k);

    // Upper bound of the index
    T c = ceil(k);

    // If the upper and lower bounds are equal,
    // we can calculate the percentile directly
    // by the index k
    if(fabs(f - c) < datadefs::EPS) {
      prc = x[static_cast<size_t>(k)];
    } else {

      // Otherwise we will interpolate linearly based on the
      // distances from the intermediate point (k) to both
      // bounds: ceil->k and k->floor
      T d0 = x[static_cast<size_t>(f)] * (c - k);
      T d1 = x[static_cast<size_t>(c)] * (k - f);

      // This operation equals to the weighted average,
      // which in other words is the interpolated percentile
      // we were after
      prc = d0 + d1;
    }

    // Finally return the calculated percentile
    return( prc );
    
  }

  void transformLogistic(size_t nCategories, vector<num_t>& prediction, vector<num_t>& probability);

  /**
     Error function
     NOTE: see http://en.wikipedia.org/wiki/Error_function
  */
  num_t erf(num_t x);

  void adjustPValues(vector<num_t>& pValues, const size_t nTests);

  /**
     Two-sample t-test
     NOTE: see http://en.wikipedia.org/wiki/Student's_t-test
  */
  num_t ttest(const vector<num_t>& x,
              const vector<num_t>& y,
	      const bool WS = false);

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
  unordered_map<T,size_t> frequency(const vector<T>& x) {
    unordered_map<T,size_t> freq;
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
    unordered_map<T,size_t> freq = frequency(x);
    typename unordered_map<T,size_t>::const_iterator maxElement( freq.begin() );
    for ( typename unordered_map<T,size_t>::const_iterator it(freq.begin()); it != freq.end(); ++it ) {
      if ( it->second > maxElement->second ) {
	maxElement = it;
      }
    }
    return( maxElement->first );
  }

  template<typename T>
  size_t nMismatches(const vector<T>& x, const T& y) {
    size_t count = 0;
    for ( size_t i = 0; i < x.size(); ++i ) {
      if ( x[i] != y ) {
	++count;
      }
    }
    return( count );
  }
  
  template<typename T>
  map<T,map<T,size_t> > confusionMap(const vector<T>& x, const vector<T>& y) {

    assert(x.size() == y.size());

    map<T,map<T,size_t> > cMap;

    set<T> allClasses;

    for ( size_t i = 0; i < x.size(); ++i ) {
      T a = x[i];
      T b = y[i];
      allClasses.insert(a);
      allClasses.insert(b);
      if ( cMap[a].find(b) == cMap[a].end() ) {
	cMap[a][b] = 1;
      } else {
	++cMap[a][b];
      }
    }

  }

  template<typename T> 
  num_t categoricalError(const vector<T>& x, const vector<T>& y) {

    assert( x.size() == y.size() );

    size_t n = x.size();

    if ( n == 0 ) {
      return(datadefs::NUM_NAN);
    }

    num_t ret = 0.0;

    for ( size_t i = 0; i < n; ++i ) {
      ret += static_cast<num_t>( x[i] != y[i] ) / n;
    }
    
    return( ret );
    
  }
  
  num_t numericalError(const vector<num_t>& x, const vector<num_t>& y);


  num_t gamma(const vector<num_t>& x, const size_t nCategories); 
  
  //num_t squaredError(const vector<num_t>& x);
  
  //num_t squaredError(const vector<num_t>& x, const num_t mu); 

  // Unbiased variance estimate: 1/(n-1)*sum(y-y_i)^2
  num_t var(const vector<num_t>& x);

  num_t var(const vector<num_t>& x, const num_t& mu);

  /**
     Updates the squared frequency by ADDING x_n to the set
     NOTE: NANs will corrupt the data
  */
  template<typename T>
  inline void incrementSquaredFrequency(const T& x_n,
					unordered_map<T,size_t>& freq,
					size_t& sqFreq) {


    // Check if the value already exists in the frequency map
    typename unordered_map<T,size_t>::iterator it(freq.find(x_n));
    if(it == freq.end()) {

      // If not, squared frequency becomes updated by 1
      sqFreq += 1;
      freq[x_n] = 1;

    } else {

      // Otherwise the squared frequency becomes updated by 
      // 2*freq + 1
      sqFreq += 2*freq[x_n] + 1;
      it->second++; //freq[x_n];

    }
  }
 
  /**
     Updates the squared frequency by REMOVING x_n
     from the set
     NOTE: NANs will corrupt the data
  */ 
  template<typename T>
  inline void decrementSquaredFrequency(const T& x_n,
                                        unordered_map<T,size_t>& freq,
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

    //cout << - 1.0 * sf_tot   / ( 1.0 * n_tot * n_tot   ) << " + " << 1.0 * sf_left  / ( 1.0 * n_tot * n_left  ) << " + " << 1.0 * sf_right / ( 1.0 * n_tot * n_right ) << endl;

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
