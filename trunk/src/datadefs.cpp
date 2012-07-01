#include "datadefs.hpp"
#include<math.h>
#include<cassert>
#include<sstream>
#include<algorithm>
#include<iostream>
#include<limits>

#include "datadefs.hpp"
#include "utils.hpp" // This will be removed after all utilities currently under datadefs are properly relocated

using namespace std;

////////////////////////////////////////////////////////////
// CONSTANTS
////////////////////////////////////////////////////////////
const datadefs::num_t datadefs::NUM_NAN = numeric_limits<double>::quiet_NaN();//numeric_limits<double>::infinity();
const string datadefs::STR_NAN = "NA";
const datadefs::num_t datadefs::NUM_INF = numeric_limits<double>::infinity();
const size_t datadefs::MAX_IDX = numeric_limits<size_t>::max();
const datadefs::num_t datadefs::EPS = 1e-18; //1e-12;
const datadefs::num_t datadefs::PI = 3.1415926535;
const datadefs::num_t datadefs::A = 0.140012;
const datadefs::num_t datadefs::LOG_OF_MAX_NUM = 70.0; /** !! Potentially
                                                        * spurious. Do you mean
                                                        * the log of the
                                                        * maximum number
                                                        * expressible as a
                                                        * num_t? */

// List of NaN's adapted from http://en.wikipedia.org/wiki/NaN#Display
const string initNANs[] = {"NA","NAN","NAN%","NANQ","NANS","QNAN","SNAN","1.#SNAN","1.#QNAN","-1.#IND","NULL","?"}; 

const set<datadefs::NAN_t> datadefs::NANs(initNANs,initNANs+12);

const string datadefs::CONTRAST = "CONTRAST";

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


////////////////////////////////////////////////////////////
// METHOD IMPLEMENTATIONS
////////////////////////////////////////////////////////////

/**
 * Convert a string vector into a category vector
 */
void datadefs::strv2catv(const vector<string>& strvec,
                         vector<datadefs::num_t>& catvec, 
			 map<string,num_t>& mapping,
			 map<num_t,string>& backMapping) {
  
  size_t n = strvec.size();
  catvec.resize(n);

  mapping.clear();
  backMapping.clear();

  num_t val = 0.0;

  //Map unique strings to values and store values in catvec as doubles 
  for(size_t strIdx = 0; strIdx < n; ++strIdx) {

    //If the string is not NaN ...
    if(!datadefs::isNAN_STR(strvec[strIdx])) {
      map<string,num_t>::iterator it;

      //Try to find the string in the map. If it's not found, extend the map...
      it = mapping.find(strvec[strIdx]);
      if(it == mapping.end()) {
        mapping.insert(pair<string,num_t>(strvec[strIdx],val));
	backMapping.insert(pair<num_t,string>(val,strvec[strIdx]));
	catvec[strIdx] = val;
        val += 1.0;
      } else {
	catvec[strIdx] = it->second; 
      }

    } else {    //If the string is defined to NaN, however...
      catvec[strIdx] = datadefs::NUM_NAN;
    }
  }  

  

  if(false) {
    cout << "mapping:" << endl;
    //for(size_t i = 0; i < strvec.size(); ++i) {
    for(map<string,num_t>::const_iterator it(mapping.begin()); it != mapping.end(); ++it) {
      cout << "mapping[" << it->first << "] => " << it->second << " => " << backMapping[ it->second ] << endl;  		    
    }
    cout << "backMapping:" << endl;
    for ( map<num_t,string>::const_iterator it(backMapping.begin()); it != backMapping.end(); ++it ) {
      cout << "backMapping[" << it->first << "] => " << it->second << " => " << mapping[ it->second ] << endl;
    }
  }

}

/**
 * Convert a string vector into a number vector
 */
void datadefs::strv2numv(const vector<string>& strvec,
                         vector<datadefs::num_t>& numvec) {
  size_t n = strvec.size();
  numvec.resize(n);

  for(size_t strIdx = 0; strIdx < n; ++strIdx) {
    //cout << strIdx << ": \"" << strvec[strIdx] << "\"" << endl;
    if(!datadefs::isNAN_STR(strvec[strIdx])) {
      numvec[strIdx] = utils::str2<datadefs::num_t>(strvec[strIdx]);
    } else {
      numvec[strIdx] = datadefs::NUM_NAN;
    }
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
 * Convert a string to a number using the intuitive solution, in base 10.

 !! Correctness: see the considerations raised here:
 http://stackoverflow.com/questions/194465/how-to-parse-a-string-to-an-int-in-c/6154614#6154614

 Refactor this method to conform to a more correct model for error handling.
*/
/*
  datadefs::num_t datadefs::str2num(const string& str) {
  
  stringstream ss( utils::chomp(str) );
  datadefs::num_t ret;
  ss >> ret;
  
  utils::streams::checkEnd(ss);
  
  return( ret );
  }
*/
  
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
 * Perform in-place sorting of the data, splitting the paired vector into
 * separate representations.
 !! Correctness: consider validations or more complete contracts for the
 methods invoked here

 !! TODO: Throw an explicit error when NaNs are passed for sorting
*/
void datadefs::sortDataAndMakeRef(const bool isIncreasingOrder,
                                  vector<num_t>& data,
                                  vector<size_t>& refIcs) {

  //assert(v.size() == ref_ics.size());
  vector<pair<num_t,size_t> > pairedv(data.size()); // !! Understandibility:
                                                    // !! consider a typedef
                                                    // !! pairedv, leaving the
                                                    // !! actual variable name
                                                    // !! as something more
                                                    // !! descriptive.

  refIcs = utils::range(data.size());

  datadefs::make_pairedv<num_t,size_t>(data,refIcs,pairedv);

  //pairedv.erase(remove_if(pairedv.begin(),pairedv.end(),&datadefs::pairedIsNAN), pairedv.end());

  if(isIncreasingOrder) {
    sort(pairedv.begin(),pairedv.end(),datadefs::increasingOrder<size_t>());
  } else {
    sort(pairedv.begin(),pairedv.end(),datadefs::decreasingOrder<size_t>());
  }

  datadefs::separate_pairedv<num_t,size_t>(pairedv,data,refIcs);
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
