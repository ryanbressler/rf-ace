#ifndef DATADEFS_HPP
#define DATADEFS_HPP

#include <cstdlib>
#include <vector>
#include <set>
#include <string>
#include <math.h>
#include <map>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include "errno.hpp"

using namespace std;

namespace datadefs {

  ////////////////////////////////////////////////////////////
  // CONSTANTS
  ////////////////////////////////////////////////////////////
  // Numerical data type
  typedef double num_t; /** Baseline numeric representation used throughout
                         *   RF-ACE. Currently, double. */

  typedef unordered_map<size_t,unordered_map<size_t,size_t> > ftable_t;

  extern const num_t NUM_NAN;       /** Numeric representation of not-a-number */
  extern const string STR_NAN;
  extern const num_t EPS;           /** Desired relative error. Literally,
                                     *   "machine EPSilon." See:
                                     *   http://en.wikipedia.org/wiki/Machine_epsilon
                                     *   and http://www.mathworks.com/help/techdoc/ref/eps.html
                                     *   for discussion
                                     */
  
  extern const num_t NUM_INF;       /** Numeric representation of positive infinity */

  extern const size_t MAX_IDX;
  extern const size_t MAX_THREADS;

  extern const num_t A;             /** Numeric constant used to estimate the
                                     *   error function of a normal distribution,
                                     *   with properties given here:
                                     *   http://homepages.physik.uni-muenchen.de/~Winitzki/erf-approx.pdf
                                     */
  extern const num_t NUM_PI;            /** Numeric representation of PI */
  extern const num_t LOG_OF_MAX_NUM;/** Numeric representation of the log of
                                     *   MAX_NUM !! ?? */

  //NaNs supported in the delimited file 
  typedef string NAN_t;             /** Used to represent NANs as a string */
  extern const set<NAN_t> NANs;     /** The complete set of string
                                     *   representations of NAN */

  extern const char tokenDelimiters[];

  extern const string CONTRAST;

  typedef enum {RF,GBT,CART,UNKNOWN} forest_t;
  extern const map<string,forest_t> forestTypeAssign;

  extern const bool       SF_DEFAULT_NO_NA_BRANCHING;

  // Random Forest default configuration
  extern const size_t     RF_DEFAULT_N_TREES;
  extern const size_t     RF_DEFAULT_M_TRY;
  extern const size_t     RF_DEFAULT_N_MAX_LEAVES;
  extern const size_t     RF_DEFAULT_NODE_SIZE;
  extern const num_t      RF_DEFAULT_IN_BOX_FRACTION;
  extern const num_t      RF_DEFAULT_SAMPLE_WITH_REPLACEMENT;
  extern const bool       RF_DEFAULT_USE_CONTRASTS;
  extern const num_t      RF_DEFAULT_CONTRAST_FRACTION;
  extern const bool       RF_DEFAULT_IS_RANDOM_SPLIT;
  extern const num_t      RF_DEFAULT_SHRINKAGE;

  // Gradient Boosting Trees default configuration
  extern const size_t     GBT_DEFAULT_N_TREES;
  extern const size_t     GBT_DEFAULT_M_TRY;
  extern const size_t     GBT_DEFAULT_N_MAX_LEAVES;
  extern const size_t     GBT_DEFAULT_NODE_SIZE;
  extern const num_t      GBT_DEFAULT_IN_BOX_FRACTION;
  extern const num_t      GBT_DEFAULT_SAMPLE_WITH_REPLACEMENT;
  extern const bool       GBT_DEFAULT_USE_CONTRASTS;
  extern const num_t      GBT_DEFAULT_CONTRAST_FRACTION;
  extern const bool       GBT_DEFAULT_IS_RANDOM_SPLIT;
  extern const num_t      GBT_DEFAULT_SHRINKAGE;

  // CART default configuration
  extern const size_t     CART_DEFAULT_N_TREES;
  extern const size_t     CART_DEFAULT_M_TRY;
  extern const size_t     CART_DEFAULT_N_MAX_LEAVES;
  extern const size_t     CART_DEFAULT_NODE_SIZE;
  extern const num_t      CART_DEFAULT_IN_BOX_FRACTION;
  extern const num_t      CART_DEFAULT_SAMPLE_WITH_REPLACEMENT;
  extern const bool       CART_DEFAULT_USE_CONTRASTS;
  extern const num_t      CART_DEFAULT_CONTRAST_FRACTION;
  extern const bool       CART_DEFAULT_IS_RANDOM_SPLIT;
  extern const num_t      CART_DEFAULT_SHRINKAGE;

  // Statistical test default configuration
  extern const size_t     FILTER_DEFAULT_N_PERMS;
  extern const num_t      FILTER_DEFAULT_P_VALUE_THRESHOLD;
  extern const bool       FILTER_DEFAULT_IS_ADJUSTED_P_VALUE;
  extern const num_t      FILTER_DEFAULT_IMPORTANCE_THRESHOLD;
  extern const bool       FILTER_NORMALIZE_IMPORTANCE_VALUES;
  extern const bool       FILTER_DEFAULT_REPORT_NONEXISTENT_FEATURES;

  // Default general configuration
  extern const bool       GENERAL_DEFAULT_PRINT_HELP;
  extern const char       GENERAL_DEFAULT_DATA_DELIMITER;
  extern const char       GENERAL_DEFAULT_HEADER_DELIMITER;
  extern const size_t     GENERAL_DEFAULT_MIN_SAMPLES;
  extern const int        GENERAL_DEFAULT_SEED;
  extern const size_t     GENERAL_DEFAULT_N_THREADS;
  extern const bool       GENERAL_DEFAULT_IS_MAX_THREADS;
  extern const num_t      GENERAL_DEFAULT_FEATURE_WEIGHT;

  
  ////////////////////////////////////////////////////////////
  // METHOD DECLARATIONS
  ////////////////////////////////////////////////////////////
  string toUpperCase(const string& str);

  bool isInteger(const string& str, int& integer);
  
  void countRealValues(vector<num_t> const& data, size_t& nRealValues);

  void map_data(vector<num_t> const& data, 
                unordered_map<num_t,vector<size_t> >& datamap,
                size_t& nRealValues);

  ////////////////////////////////////////////////////////////
  // INLINE METHOD DEFINITIONS
  ////////////////////////////////////////////////////////////
  /**
   * Checks an input string against a dictionary of known representations of
   *  not-a-number. Returns true if this string exactly contains one of these
   *  representations; false otherwise.
   */
  inline bool isNAN_STR(const string& str) {
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

}

#endif
