#include "datadefs.hpp"
#include <math.h>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <limits>

#ifndef NOTHREADS
#include <thread>
#endif


using namespace std;

////////////////////////////////////////////////////////////
// CONSTANTS
////////////////////////////////////////////////////////////
const datadefs::num_t datadefs::NUM_NAN = numeric_limits<double>::quiet_NaN();//numeric_limits<double>::infinity();
const string datadefs::STR_NAN = "NA";
const datadefs::num_t datadefs::NUM_INF = numeric_limits<double>::infinity();
const size_t datadefs::MAX_IDX = numeric_limits<int32_t>::max() - 1;
const datadefs::num_t datadefs::EPS = 1e-18; //1e-12;
const datadefs::num_t datadefs::NUM_PI = 3.1415926535;
const datadefs::num_t datadefs::A = 0.140012;
const datadefs::num_t datadefs::LOG_OF_MAX_NUM = 70.0; /** !! Potentially
                                                        * spurious. Do you mean
                                                        * the log of the
                                                        * maximum number
                                                        * expressible as a
                                                        * num_t? */

// List of NaN's adapted from http://en.wikipedia.org/wiki/NaN#Display
const string initNANs[] = {"NA","NAN","NAN%","NANQ","NANS","QNAN","SNAN","1.#SNAN","1.#QNAN","-1.#IND","NULL","?"}; 

const char datadefs::tokenDelimiters[] = " \t,.;:?!@'\"-\n";

const set<datadefs::NAN_t> datadefs::NANs(initNANs,initNANs+12);

const string datadefs::CONTRAST = "CONTRAST";

#ifndef NOTHREADS
const size_t datadefs::MAX_THREADS = thread::hardware_concurrency();
#endif

#ifdef NOTHREADS 
const size_t datadefs::MAX_THREADS = 1;
#endif

// Random Forest default configuration
const size_t          datadefs::RF_DEFAULT_N_TREES = 100;
const size_t          datadefs::RF_DEFAULT_M_TRY = 0;
const size_t          datadefs::RF_DEFAULT_N_MAX_LEAVES = datadefs::MAX_IDX;
const size_t          datadefs::RF_DEFAULT_NODE_SIZE = 3;
const datadefs::num_t datadefs::RF_DEFAULT_IN_BOX_FRACTION = 1.0;
const datadefs::num_t datadefs::RF_DEFAULT_SAMPLE_WITH_REPLACEMENT = true;
const bool            datadefs::RF_DEFAULT_USE_CONTRASTS = false;
const datadefs::num_t datadefs::RF_DEFAULT_CONTRAST_FRACTION = 0.5;
const bool            datadefs::RF_DEFAULT_IS_RANDOM_SPLIT = true;
const datadefs::num_t datadefs::RF_DEFAULT_SHRINKAGE = 0.0;

// Gradient Boosting Trees default configuration
const size_t          datadefs::GBT_DEFAULT_N_TREES = 100;
const size_t          datadefs::GBT_DEFAULT_M_TRY = 0;
const size_t          datadefs::GBT_DEFAULT_N_MAX_LEAVES = 6;
const size_t          datadefs::GBT_DEFAULT_NODE_SIZE = 3;
const datadefs::num_t datadefs::GBT_DEFAULT_IN_BOX_FRACTION = 0.5;
const datadefs::num_t datadefs::GBT_DEFAULT_SAMPLE_WITH_REPLACEMENT = false;
const bool            datadefs::GBT_DEFAULT_USE_CONTRASTS = false;
const datadefs::num_t datadefs::GBT_DEFAULT_CONTRAST_FRACTION = 0.5;
const bool            datadefs::GBT_DEFAULT_IS_RANDOM_SPLIT = false;
const datadefs::num_t datadefs::GBT_DEFAULT_SHRINKAGE = 0.1;

// CART default configuration
const size_t          datadefs::CART_DEFAULT_N_TREES = 1;
const size_t          datadefs::CART_DEFAULT_M_TRY = 0;
const size_t          datadefs::CART_DEFAULT_N_MAX_LEAVES = 0;
const size_t          datadefs::CART_DEFAULT_NODE_SIZE = 3;
const datadefs::num_t datadefs::CART_DEFAULT_IN_BOX_FRACTION = 1.0;
const datadefs::num_t datadefs::CART_DEFAULT_SAMPLE_WITH_REPLACEMENT = false;
const bool            datadefs::CART_DEFAULT_USE_CONTRASTS = false;
const datadefs::num_t datadefs::CART_DEFAULT_CONTRAST_FRACTION = 0.5;
const bool            datadefs::CART_DEFAULT_IS_RANDOM_SPLIT = false;
const datadefs::num_t datadefs::CART_DEFAULT_SHRINKAGE = 0;

// Statistical test default configuration
const size_t          datadefs::FILTER_DEFAULT_N_PERMS = 20;
const datadefs::num_t datadefs::FILTER_DEFAULT_P_VALUE_THRESHOLD = 0.05;
const bool            datadefs::FILTER_DEFAULT_IS_ADJUSTED_P_VALUE = false;
const datadefs::num_t datadefs::FILTER_DEFAULT_IMPORTANCE_THRESHOLD = 10;
const bool            datadefs::FILTER_NORMALIZE_IMPORTANCE_VALUES = false;
const bool            datadefs::FILTER_DEFAULT_REPORT_NONEXISTENT_FEATURES = false;


// Default general configuration
const bool            datadefs::GENERAL_DEFAULT_PRINT_HELP = false;
const char            datadefs::GENERAL_DEFAULT_DATA_DELIMITER = '\t';
const char            datadefs::GENERAL_DEFAULT_HEADER_DELIMITER = ':';
const size_t          datadefs::GENERAL_DEFAULT_MIN_SAMPLES = 5;
const int             datadefs::GENERAL_DEFAULT_SEED = -1;
const size_t          datadefs::GENERAL_DEFAULT_N_THREADS = 1;
const bool            datadefs::GENERAL_DEFAULT_IS_MAX_THREADS = false;
const datadefs::num_t datadefs::GENERAL_DEFAULT_FEATURE_WEIGHT = 0;


////////////////////////////////////////////////////////////
// HELPER FUNCTIONS
////////////////////////////////////////////////////////////

/**
 * Promote each character in a sequence to uppercase. Effectively, a wrapper
 * around std::transform.
 */

string datadefs::toUpperCase(const string& str) {
  int (*pf)(int) = toupper;  
  string strcopy(str);
  transform(strcopy.begin(), strcopy.end(), strcopy.begin(), pf);
  return(strcopy);
}

bool datadefs::isInteger(const string& str, int& integer) {
  stringstream ss(str);
  if(ss >> integer && ss.eof()) {
    return(true);
  } else {
    return(false);
  }
}


bool datadefs::is_unique(const vector<string>& strvec) {

  set<string> strset;
  for ( size_t i = 0; i < strvec.size(); ++i ) {
    if ( strset.find(strvec[i]) != strset.end() ) {
      return false;
    } else {
      strset.insert(strvec[i]);
    }
  }
  return true;
}

  
/**
 * Determines the cardinality of a given input data set
 !! Input sanitization: contains no checks for NaN.
*/
void datadefs::cardinality(const vector<datadefs::num_t>& data, size_t& cardinality) {
  set<datadefs::num_t> categories;
  for(size_t i = 0; i < data.size(); ++i) {
    if (!datadefs::isNAN(data[i])) { // Filter out NaN, as it causes unintended results
                              //  in std::set::insert.
      categories.insert(data[i]);
    }
  }
  cardinality = categories.size();

}

/**
 * Count all values that aren't transfinite
 !! Correctness: what about representations of infinity? And to be entirely
 pedantic: signaling NaN, post-trap? These should have specific non-guarantees.
*/
void datadefs::countRealValues(vector<num_t> const& data, size_t& nRealValues) {
  nRealValues = 0;
  for(size_t i = 0; i < data.size(); ++i) {
    if(!datadefs::isNAN(data[i])) {
      ++nRealValues;
    }
  }  
}


/**
   !! Documentation
*/
void datadefs::map_data(vector<datadefs::num_t> const& data, 
                        map<datadefs::num_t,vector<size_t> >& datamap, 
                        size_t& nRealValues) {
  datamap.clear();
  map<datadefs::num_t,vector<size_t> >::iterator it;
  nRealValues = 0;
  for(size_t i = 0; i < data.size(); ++i) {
    if(!datadefs::isNAN(data[i])) {
      ++nRealValues;
      it = datamap.find(data[i]);
      if(it == datamap.end()) {
        datamap.insert(pair<datadefs::num_t,vector<size_t> >(data[i],vector<size_t>(1,i)));
      } else {
        it->second.push_back(i);
      }
    }
  }
}
