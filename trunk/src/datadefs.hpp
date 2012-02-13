#ifndef DATADEFS_HPP
#define DATADEFS_HPP

#include<cstdlib>
#include<vector>
#include<set>
#include<string>
#include<math.h>
#include<map>
#include<cassert>
#include<iostream>
#include<algorithm>
#include "errno.hpp"

using namespace std;

namespace datadefs {

  ////////////////////////////////////////////////////////////
  // CONSTANTS
  ////////////////////////////////////////////////////////////
  // Numerical data type
  typedef double num_t; /** Baseline numeric representation used throughout
                         *   RF-ACE. Currently, double. */

  // Data type for mapping string literals to numbers
  //typedef map<string,num_t> fwMapping_t;

  // Data type for mapping numbers (back) to string literals
  //typedef map<num_t,string> bwMapping_t;

  extern const num_t NUM_NAN;       /** Numeric representation of not-a-number */
  extern const string STR_NAN;
  extern const num_t EPS;           /** Desired relative error. Literally,
                                     *   "machine EPSilon." See:
                                     *   http://en.wikipedia.org/wiki/Machine_epsilon
                                     *   and http://www.mathworks.com/help/techdoc/ref/eps.html
                                     *   for discussion
                                     */
  
  extern const num_t NUM_INF;       /** Numeric representation of positive infinity */

  extern const num_t A;             /** Numeric constant used to estimate the
                                     *   error function of a normal distribution,
                                     *   with properties given here:
                                     *   http://homepages.physik.uni-muenchen.de/~Winitzki/erf-approx.pdf
                                     */
  extern const num_t PI;            /** Numeric representation of PI */
  extern const num_t LOG_OF_MAX_NUM;/** Numeric representation of the log of
                                     *   MAX_NUM !! ?? */

  //NaNs supported in the delimited file 
  typedef string NAN_t;             /** Used to represent NANs as a string */
  extern const set<NAN_t> NANs;     /** The complete set of string
                                     *   representations of NAN */

  
  ////////////////////////////////////////////////////////////
  // METHOD DECLARATIONS
  ////////////////////////////////////////////////////////////
  string toUpperCase(const string& str);

  bool isInteger(const string& str, int& integer);
  
  void strv2catv(const vector<string>& strvec, vector<num_t>& catvec, map<string,num_t>& mapping, map<num_t,string>& backMapping);
  void strv2numv(const vector<string>& strvec, vector<num_t>& numvec);
  num_t str2num(const string& str);

  //string chomp(const string& str);

  bool is_unique(const vector<string>& strvec);

  //void mean(vector<num_t> const& data, num_t& mu, size_t& nreal);

  // Newly added function set which used in calculating leaf node predictions
  void mean(vector<num_t> const& data, num_t& mu, size_t& numClasses);
  //void mode(const vector<num_t>& data, num_t& mode, const size_t numClasses);
  //void gamma(const vector<num_t>& data, num_t& gamma, const size_t numClasses);

  void cardinality(const vector<num_t>& data, size_t& cardinality);

  void sqerr(vector<num_t> const& data, 
             num_t& mu, 
             num_t& se,
             size_t& nreal);      

  void countRealValues(vector<num_t> const& data, size_t& nRealValues);

  void count_freq(vector<num_t> const& data, map<num_t,size_t>& cat2freq, size_t& nRealValues);

  void map_data(vector<num_t> const& data, 
                map<num_t,vector<size_t> >& datamap,
                size_t& nRealValues);

  void gini(vector<num_t> const& data,
            num_t& giniIndex,
            size_t& nRealValues);

  void gini(map<num_t,size_t> const& cat2freq,
            num_t& giniIndex);

  void sqfreq(vector<num_t> const& data, 
              map<num_t,size_t>& freq, 
              size_t& sqFreq, 
              size_t& nRealValues);

  void forward_sqfreq(const num_t x_n,
                      size_t& n,
                      map<num_t,size_t>& freq,
                      size_t& sqFreq);
  
  void forward_backward_sqfreq(const num_t x_n,
                               size_t& n_left,
                               map<num_t,size_t>& freq_left,
                               size_t& sf_left,
                               size_t& n_right,
                               map<num_t,size_t>& freq_right,
                               size_t& sf_right);
  
  void range(vector<size_t>& ics);
  void sortDataAndMakeRef(const bool isIncreasingOrder, vector<num_t>& data, vector<size_t>& refIcs);

  void utest(vector<num_t> const& x,
             vector<num_t> const& y,
             num_t& pvalue);

  num_t erf(num_t x);

  void pearson_correlation(vector<num_t> const& x,
                           vector<num_t> const& y,
                           num_t& corr);

  num_t ttest(vector<num_t> const& x, 
	      vector<num_t> const& y);


  num_t regularized_betainc(const num_t x, 
			    const size_t a);


  
  ////////////////////////////////////////////////////////////
  // INLINE METHOD DEFINITIONS
  ////////////////////////////////////////////////////////////
  /**
   * Checks an input string against a dictionary of known representations of
   *  not-a-number. Returns true if this string exactly contains one of these
   *  representations; false otherwise.
   */
  inline bool isNAN(const string& str) {
    return( NANs.find(toUpperCase(str)) != NANs.end() ? true : false );
  }


  /**
   * Performs an equivalence test to discern if this value is NAN.
   */
  inline bool isNAN(const num_t value) {
    return( value != value ? true : false );
  }

  /**
   * Checks if a data array contains at least one representation of NAN
   */
  inline bool containsNAN(const vector<num_t>& data) {
    
    for(size_t i = 0; i < data.size(); ++i) {
      if(data[i] != data[i]) { return(true); }
    }
    return(false);
  }

  inline bool pairedIsNAN(const pair<num_t,size_t>& value) {
    return( value.first != value.first ? true : false );
  }

  /**
     !! Document
     * Calculates the squared error, given x, n, mu, and standard error
     */
  inline void forward_sqerr(const num_t& x_n,
                            size_t& n,
                            num_t& mu,
                            num_t& se) {  

    if(datadefs::isNAN(x_n)) { return; }  // Check for NAN
    
    ++n; if (n == 0) { throw ERRNO_NUMERIC_OVERFLOW; }
    num_t mu_old = mu;
    mu += (x_n - mu) / n;
    
    //If there are already at least two data points, squared error can be calculated, otherwise assign se_left := 0.0
    if(n > 1) {
      se += (x_n - mu) * (x_n - mu_old);
    } else {
      
      se = 0.0; /** Implementation note: this may spuriously invoke on
                     overflow. size_t is assumed to always be unsigned in
                     accordance with the 1999 ISO C standard (C99). */
    }
  } 

  /**
     !! Document
     * Calculates the squared error from both supplied branches, x and each
     *  branch's respective n, mu, and standard error.
     */
  inline void forward_backward_sqerr(const num_t& x_n,
                                     size_t& n_left,
                                     num_t& mu_left,
                                     num_t& se_left,
                                     size_t& n_right,
                                     num_t& mu_right,
                                     num_t& se_right) {
    
    if(datadefs::isNAN(x_n)) { return; }  // Check for NAN
    
    //assert(n_right > 0);
    ++n_left;   if (n_left == 0)                        { throw ERRNO_NUMERIC_OVERFLOW; }
    --n_right;  if (n_right == static_cast<size_t>(-1)) { throw ERRNO_NUMERIC_UNDERFLOW; }
    
    // As long as there are at least two data points on the "right" branch,
    //  squared error can be calculated, otherwise assign se_right := 0.0
    
    if(n_right > 1) {
      datadefs::num_t mu_old = mu_right;
      mu_right -= (x_n - mu_right) / n_right;
      se_right -= (x_n - mu_right) * (x_n - mu_old);
    } else if(n_right == 1) {
      mu_right -= (x_n - mu_right) / n_right;
      se_right = 0.0;
    } else {
      mu_right = 0.0;
      se_right = 0.0;
    }
    
    // Add x_n to "left" and update mean and squared error
    datadefs::num_t mu_old = mu_left;
    mu_left += (x_n - mu_left) / n_left;
    
    // If there are already at least two data points on the "left" branch,
    //  squared error can be calculated, otherwise assign se_left := 0.0
    
    if(n_left > 1) { se_left += (x_n - mu_left) * (x_n - mu_old); } else { se_left = 0.0; }
  }

  
  /**
   * A comparator functor that can be passed to STL::sort. Assumes that one is
   *  comparing first elements of pairs, first type being num_t and second T.
   */
  template <typename T> struct increasingOrder {
    bool operator ()(pair<datadefs::num_t,T> const& a, pair<datadefs::num_t,T> const& b) {
      return(a.first < b.first);
    }
  }; 

  /**
   * A comparator functor that can be passed to STL::sort. Assumes that one is
   *  comparing first elements of pairs, first type being num_t and second T.
   */
  template <typename T> struct decreasingOrder {
    bool operator ()(
      pair<datadefs::num_t,T> const& a,
      pair<datadefs::num_t,T> const& b
      ) {
      return(a.first > b.first);
    }
  };

  /**
   * A comparator functor that can be passed to STL::sort. Assumes that one is
   *  comparing second elements of maps, first type being num_t and second size_t.
   */
  struct freqIncreasingOrder {
    bool operator ()(
      const map<datadefs::num_t,size_t>::value_type& a,
      const map<datadefs::num_t,size_t>::value_type& b
      ) {
      return(a.second < b.second);
    }
  };

  /**
     !! Document
     * Union and splitting operations for vectors. Consider fully documenting.
     !! Correctness: This may cause problems if T1 or T2 are declared as unsafe
         for assignment. 
     */
  template <typename T1,typename T2> void make_pairedv(vector<T1> const& v1,
                                                       vector<T2> const& v2,
                                                       vector<pair<T1,T2> >& p) {
    size_t n = v1.size();
    assert(n == v2.size());
    p.resize(n);
    for(size_t i = 0; i < n; ++i) {
      p[i].first = v1[i]; p[i].second = v2[i];
    }
  }

  
  /**
     !! Document
     * Union and splitting operations for vectors. Consider documenting.
     */
  template <typename T1,typename T2> void separate_pairedv(vector<pair<T1,T2> > const& p,
                                                           vector<T1>& v1,
                                                           vector<T2>& v2) {
    //assert(v1.size() == v2.size() && v2.size() == p.size());
    size_t n = p.size();
    v1.resize(n);
    v2.resize(n);
    for(size_t i = 0; i < n; ++i) {
      v1[i] = p[i].first;
      v2[i] = p[i].second;
    }
  }

  /**
   * Sorts a given input data vector of type T based on a given reference
   * ordering of type vector<int>.
   !! Correctness: this will fail if any of the contents of refIcs fall outside
       of the normal scope of vector<T>& data.
   */
  template <typename T> void sortFromRef(vector<T>& data,
					 vector<size_t> const& refIcs
					 ) {
    assert(data.size() == refIcs.size());  
    vector<T> foo = data;
    int n = data.size();
    for (int i = 0; i < n; ++i) {
      data[i] = foo[refIcs[i]];
    }
  }

  /**
   * Prints a vector
   */
  template <typename T> void print(const vector<T>& x) {
    for(size_t i = 0; i < x.size(); ++i) {
      cout << " " << x[i];
    }
    cout << endl;
  }

  template <typename T> void print(const set<T>& x) {
    
    for( typename set<T>::const_iterator it = x.begin(); it != x.end(); ++it) {
      cout << " " << *it;
    }
    cout << endl;
  }

}

#endif
