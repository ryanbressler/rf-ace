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
const datadefs::num_t datadefs::EPS = 1e-18; //1e-12;
const datadefs::num_t datadefs::PI = 3.1415926535;
const datadefs::num_t datadefs::A = 0.140012;
const datadefs::num_t datadefs::LOG_OF_MAX_NUM = 70.0; /** !! Potentially
                                                        * spurious. Do you mean
                                                        * the log of the
                                                        * maximum number
                                                        * expressible as a
                                                        * num_t? */

const string initNANs[] = {"NA","NAN","?"}; /** !! Incomplete definition:
                                             * consider including:
                                             *  http://en.wikipedia.org/wiki/NaN#Display */
const set<datadefs::NAN_t> datadefs::NANs(initNANs,initNANs+3);


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
      numvec[strIdx] = str2num(strvec[strIdx]);
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
datadefs::num_t datadefs::str2num(const string& str) {

  // Initialize and use the stringstream
  //stringstream ss(terminatorIdx != -1 ? str.substr(0,terminatorIdx) : str);
  stringstream ss( utils::chomp(str) );
  datadefs::num_t ret;
  ss >> ret;
  
  if (ss.fail()) {  // Ensure reading didn't fail
    cerr << "datadefs::str2num: ERROR: parameter '" << str
	 << "' could not be read properly. Quitting... "
	 << endl;
    assert(false);
    ret = 0.0;
    
  } else if (!ss.eof()) {   // Ensure eofbit is set
    cerr << "datadefs::str2num: ERROR: parameter '" << str
	 << "' was only partially read. Quitting... "
	 << endl;
    assert(false);
  }
  return(ret);
}

/**
 * Takes the arithmetic mean of a vector of input values
 !! Side effects: mutates nRealValues and mu regardless of the final
 result. This applies to most methods defined in RF-ACE, actually.

 !! Error sieve: silently ignores NaN

 !! Consistency: consider replacing "mu" with "mean" in the function signature
     to make this consistent with the mode and gamma values. Yes, mu is the
     common first-order notation for the mean, but explicitly calling it out
     is more expressive.

 !! Consistency: what's the meaningful difference (pun!) between this and the
     now-since-commented-out datadefs::mean? This reads like a function derived
     from shotgun debugging.


*/
void datadefs::mean(vector<datadefs::num_t> const& data, datadefs::num_t& mu, size_t& nRealValues) {
 
  mu = 0.0;
  nRealValues = 0;
 
  for(size_t i = 0; i < data.size(); ++i) {
    if(!datadefs::isNAN(data[i])) {
      ++nRealValues;
      mu += data[i];
    }
  }

  if(nRealValues > 0) {
    mu /= nRealValues;
  }
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
 * Calculates the squared error
 !! Possible duplication: how does this differ with similar inline methods in datadefs.hpp?
*/
/*
  void datadefs::sqerr(vector<datadefs::num_t> const& data, 
  datadefs::num_t& mu, 
  datadefs::num_t& se,
  size_t& nRealValues) {
  
  datadefs::mean(data,mu,nRealValues);
  
  se = 0.0;
  for(size_t i = 0; i < data.size(); ++i) {
  if(!datadefs::isNAN(data[i])) {
  se += pow(data[i] - mu,2);
  }
  }
  }
*/

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
  refIcs.resize(data.size()); // The actual allocation of refIcs is irrelevant;
                              //  its contained data will be flattened and
                              //  resized. Note that any values that are unsafe
                              //  for equals may cause an unintended memory leak.
  datadefs::range(refIcs);

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
 * Gini coefficient (see: http://en.wikipedia.org/wiki/Gini_coefficient)
 !! Documentation
*/
/*
  void datadefs::gini(vector<datadefs::num_t> const& data,
  datadefs::num_t& giniIndex,
  size_t& nRealValues) {
  map<datadefs::num_t,size_t> freq;
  datadefs::count_freq(data,freq,nRealValues);
  datadefs::gini(freq,giniIndex);
  }
*/

/**
 * Gini coefficient (see: http://en.wikipedia.org/wiki/Gini_coefficient)
 !! Documentation
*/
/*
  void datadefs::gini(map<datadefs::num_t,size_t> const& cat2freq, 
  datadefs::num_t& giniIndex) {
  giniIndex = 0.0;
  size_t n = 0;
  map<datadefs::num_t,size_t>::const_iterator it;
  for(it = cat2freq.begin(); it != cat2freq.end(); ++it) {
  size_t freq_new = it->second;
  giniIndex += freq_new * freq_new;
  n += freq_new;
  }
  if(n) {
  giniIndex = 1 - giniIndex / ( n*n );
  }
  }
*/


/**
   !! Documentation
   !! Critical mutator of the cat2freq map. This should be renamed; count_freq
   !! sounds like an innocuous accessor that counts the elements in a given
   !! collection of descriptive type "freq".
*/
/*
  void datadefs::count_freq(vector<datadefs::num_t> const& data, map<datadefs::num_t,size_t>& cat2freq, size_t& nRealValues) {
  cat2freq.clear();
  map<datadefs::num_t,size_t>::const_iterator it;
  nRealValues = 0;
  for(size_t i = 0; i < data.size(); ++i) {
  if(!datadefs::isNAN(data[i])) {
  ++nRealValues;
  it = cat2freq.find(data[i]);
  if(it == cat2freq.end()) {
  cat2freq.insert(pair<datadefs::num_t,size_t>(data[i],1));
  } else {
  ++cat2freq[data[i]];
  }
  }
  }
  }
*/

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
        vector<size_t> foo;
        foo.push_back(i);
        datamap.insert(pair<datadefs::num_t,vector<size_t> >(data[i],foo));
      } else {
        it->second.push_back(i);
      }
    }
  }
}

// !! Documentation: this is just a drop-in replacement for Python's range()
// !! function, hinted by the size of the input vector. It mutates ics,
// !! specifying a 0-based range in all cases, and could be made more robust if
// !! the starting value could be given.
void datadefs::range(vector<size_t>& ics) {
  for(size_t i = 0; i < ics.size(); ++i) {
    ics[i] = i;
  }
}
