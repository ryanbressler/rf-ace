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
//using datadefs::ForestType;

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

const char datadefs::tokenDelimiters[] = " \t,.;:?!@'\"-\n\0";

// List of NaN's adapted from http://en.wikipedia.org/wiki/NaN#Display
const set<datadefs::NAN_t> datadefs::NANs = {"NA","NAN","NULL","?"};

const string datadefs::CONTRAST = "CONTRAST";

#ifndef NOTHREADS
const size_t datadefs::MAX_THREADS = thread::hardware_concurrency();
#endif

#ifdef NOTHREADS 
const size_t datadefs::MAX_THREADS = 1;
#endif

//enum ForestType {RF, GBT, CART, UNKNOWN};

const map<string,datadefs::forest_t> datadefs::forestTypeAssign = { {"RF",datadefs::forest_t::RF}, {"GBT",datadefs::forest_t::GBT}, {"CART",datadefs::forest_t::CART} };

const bool                      datadefs::SF_DEFAULT_NO_NA_BRANCHING = false;
const vector<datadefs::num_t>   datadefs::SF_DEFAULT_QUANTILES = {};

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
const size_t          datadefs::CART_DEFAULT_N_MAX_LEAVES = datadefs::MAX_IDX;
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
                        unordered_map<datadefs::num_t,vector<size_t> >& datamap, 
                        size_t& nRealValues) {

  datamap.clear();
  datamap.reserve(data.size());

  unordered_map<datadefs::num_t,vector<size_t> >::iterator it;
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
