#include "datadefs.hpp"
#include<math.h>
#include<cassert>
#include<sstream>
#include<algorithm>
#include<iostream>
#include<limits>

using namespace std;

////////////////////////////////////////////////////////////
// CONSTANTS
////////////////////////////////////////////////////////////
const datadefs::num_t datadefs::NUM_NAN = numeric_limits<double>::quiet_NaN();//numeric_limits<double>::infinity();
const string datadefs::STR_NAN = "NA";
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


////////////////////////////////////////////////////////////
// METHOD IMPLEMENTATIONS
////////////////////////////////////////////////////////////

/**
 * Convert a string vector into a category vector
 */
void datadefs::strv2catv(vector<string>& strvec,
                         vector<datadefs::num_t>& catvec, 
			 map<string,num_t>& mapping,
			 map<num_t,string>& backMapping) {
  assert(strvec.size() == catvec.size());

  //map<string,size_t> str2valmap;

  //Reset the map
  //map<string,size_t> foo;
  //str2valmap = foo;
  mapping.clear();
  backMapping.clear();

  //mapping.insert(pair<string,num_t>("NA",datadefs::NUM_NAN));
  //backMapping.insert(pair<num_t,string>(datadefs::NUM_NAN,"NA"));

  num_t val = 0.0;

  //Map unique strings to values and store values in catvec as doubles 
  for(size_t strIdx = 0; strIdx < strvec.size(); ++strIdx) {
    //Transform string to uppercase
    //toupper(strvec[strIdx]);

    //If the string is not defined to NaN
    if(!datadefs::isNAN(strvec[strIdx])) {
      map<string,num_t>::iterator it;

      //Try to find the string in the map. If it's not found, add the map...
      it = mapping.find(strvec[strIdx]);
      if(it == mapping.end()) {
        mapping.insert(pair<string,num_t>(strvec[strIdx],val));
	backMapping.insert(pair<num_t,string>(val,strvec[strIdx]));
	catvec[strIdx] = val;
        val += 1.0;
      } else {
	catvec[strIdx] = it->second; //mapping[ strvec[strIdx] ];
      }
      //...and use the map to set the value for the output vector
      //catvec[strIdx] = val; //backMapping[strvec[strIdx]];
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
void datadefs::strv2numv(vector<string>& strvec,
                         vector<datadefs::num_t>& numvec) {
  assert(strvec.size() == numvec.size());
  
  for(size_t strIdx = 0; strIdx < strvec.size(); ++strIdx) {
    //toupper(strvec[strIdx]);
    if(!datadefs::isNAN(strvec[strIdx])) {
      numvec[strIdx] = str2num(strvec[strIdx]);
    } else {
      numvec[strIdx] = datadefs::NUM_NAN;
    }
  }
}

/**
 * Convert a string to a number using the intuitive solution, in base 10.

 !! Correctness: see the considerations raised here:
 http://stackoverflow.com/questions/194465/how-to-parse-a-string-to-an-int-in-c/6154614#6154614

 Refactor this method to conform to a more correct model for error handling.
*/
datadefs::num_t datadefs::str2num(string& str) {

  // Chop at the first newline character, if it exists
  int crIdx = str.find("\r");
  int lfIdx = str.find("\n");
  int terminatorIdx = crIdx;
  if (lfIdx != -1 && lfIdx < crIdx) {
    terminatorIdx = lfIdx;
  }

  // Initialize and use the stringstream
  stringstream ss(terminatorIdx != -1 ? str.substr(0,terminatorIdx) : str);
  datadefs::num_t ret;
  ss >> ret;
  
  if (ss.fail()) {  // Ensure reading didn't fail
    cerr << "WARNING: parameter '" << str
         << "' could not be read properly. *THIS MAY CAUSE SPURIOUS RESULTS!*"
         << endl;
    ret = 0.0;
    
  } else if (!ss.eof()) {   // Ensure eofbit is set
    cerr << "WARNING: parameter '" << str
         << "' was only partially read. *THIS MAY CAUSE SPURIOUS RESULTS!*"
         << endl;
    
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
  //cout << "sort_and_make_ref: in the beginning" << endl;
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
  //cout << "sort_and_make_ref: used range()" << endl;
  datadefs::make_pairedv<num_t,size_t>(data,refIcs,pairedv);
  //cout << "sort_and_make_ref: made pairedv" << endl;
  if(isIncreasingOrder) {
    sort(pairedv.begin(),pairedv.end(),datadefs::increasingOrder<num_t>());
  } else {
    sort(pairedv.begin(),pairedv.end(),datadefs::decreasingOrder<num_t>());
  }
  //cout << "sort_and_make_ref: pairedv sorted" << endl;
  datadefs::separate_pairedv<num_t,size_t>(pairedv,data,refIcs);
}

/**
 * Gini coefficient (see: http://en.wikipedia.org/wiki/Gini_coefficient)
 !! Documentation
*/
void datadefs::gini(vector<datadefs::num_t> const& data,
                    datadefs::num_t& giniIndex,
                    size_t& nRealValues) {
  map<datadefs::num_t,size_t> freq;
  datadefs::count_freq(data,freq,nRealValues);
  datadefs::gini(freq,giniIndex);
}

/**
 * Gini coefficient (see: http://en.wikipedia.org/wiki/Gini_coefficient)
 !! Documentation
*/
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


/**
   !! Documentation
   !! Critical mutator of the cat2freq map. This should be renamed; count_freq
   !! sounds like an innocuous accessor that counts the elements in a given
   !! collection of descriptive type "freq".
*/
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

// !! Documentation: summation function featuring the squared representation of
// !! each frequency. This should be documented.
void datadefs::sqfreq(vector<datadefs::num_t> const& data, 
                      map<datadefs::num_t,size_t>& freq, 
                      size_t& sqFreq, 
                      size_t& nRealValues) {
  sqFreq = 0;
  datadefs::count_freq(data,freq,nRealValues);
  for(map<datadefs::num_t,size_t>::const_iterator it(freq.begin()); it != freq.end(); ++it) {
    size_t freq_new = it->second;
    sqFreq += freq_new * freq_new;
  }
  // !! Correctness: this doesn't actually return anything. Instead, it mutates
  // !! sqFreq by reference, without returning any additional information.
}

// !! Correctness: mutates freq and sqFreq, does not supply an explicit
// !! contract for these values.

// !! Correctness: this isn't actually a squaring of the frequencies! This
// !! should probably be named better.

// !! Documentation: I have no idea what this is supposedly doing.
void datadefs::forward_sqfreq(const datadefs::num_t x_n,
                              size_t& n,
                              map<datadefs::num_t,size_t>& freq,
                              size_t& sqFreq) {

  if(datadefs::isNAN(x_n)) { // !! Correctness: implicit check for NAN
    return;
  }

  ++n; if (n == 0) { throw ERRNO_NUMERIC_OVERFLOW; }
  
  //cout << "sf_old=" << sf;
  //Check if the value already exists
  map<datadefs::num_t,size_t>::const_iterator it(freq.find(x_n));
  if(it == freq.end()) {
    //cout << "sf_old=" << sf;
    sqFreq += 1;
    //cout << "  sf_new=" << sf << endl;

    //If not, add a new category and set its frequency to 1...
    freq.insert(pair<datadefs::num_t,size_t>(x_n,1));

  } else {
    sqFreq += 2*freq[x_n] + 1;
    ++freq[x_n];
  }
  //cout << "  sf_new=" << sf << endl;

}

// !! Correctness: mutates a whole lot of its inputs, does not supply any
// !! explicit contract for this behavior.

// !! Correctness: this isn't actually a squaring of the frequencies! This
// !! should probably be named better.

// !! Documentation: I have no idea what this is supposedly doing.
void datadefs::forward_backward_sqfreq(const datadefs::num_t x_n,
                                       size_t& n_left,
                                       map<datadefs::num_t,size_t>& freq_left, 
                                       size_t& sf_left,
                                       size_t& n_right,
                                       map<datadefs::num_t,size_t>& freq_right,
                                       size_t& sf_right) {
 
  if(datadefs::isNAN(x_n)) {
    return;
  }

  //assert(n_right > 0);
  ++n_left;   if (n_left == 0)                        { throw ERRNO_NUMERIC_OVERFLOW; }
  --n_right;  if (n_right == static_cast<size_t>(-1)) { throw ERRNO_NUMERIC_UNDERFLOW; }

  //Check if the value already exists on left
  map<datadefs::num_t,size_t>::const_iterator it(freq_left.find(x_n));
  if(it == freq_left.end()) {
    sf_left += 1;

    //If not, add a new category and set its frequency to 1...
    freq_left.insert(pair<datadefs::num_t,size_t>(x_n,1));
      
  } else {
    sf_left += 2*freq_left[x_n] + 1;
    ++freq_left[x_n];
  }

  it = freq_right.find(x_n);
  assert(it != freq_right.end());
  assert(freq_right[x_n] > 0);

  sf_right -= 2*freq_right[x_n] - 1;
  --freq_right[x_n];

  if(freq_right[x_n] == 0) {
    freq_right.erase(x_n);
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

//DEPRECATED

// !! Documentation: performs a T-test on the input. It should be readily
// !! apparent why this exists,

// !! Spurious deprecation: this isn't used in the existing codebase, but it
// !! has value within the API. Can we consider it "deprecated?"

void datadefs::ttest(vector<datadefs::num_t> const& x, 
                     vector<datadefs::num_t> const& y, 
                     datadefs::num_t& pvalue) {
  // Sample mean and variance of x
  datadefs::num_t mean_x = 0;
  datadefs::num_t var_x = 0;
  size_t nreal_x = 0;

  datadefs::sqerr(x, mean_x, var_x, nreal_x);

  assert(nreal_x > 1);

  var_x /= (nreal_x - 1);

  // Sample mean and variance of y
  datadefs::num_t mean_y = 0;
  datadefs::num_t var_y = 0;
  size_t nreal_y = 0;

  datadefs::sqerr(y, mean_y, var_y, nreal_y);
    
  assert(nreal_y > 1);

  var_y /= (nreal_y - 1);
    
  assert(nreal_x == nreal_y);

  if(var_x < datadefs::EPS && var_y < datadefs::EPS) {
    if(fabs(mean_x - mean_y) < datadefs::EPS) {
      pvalue = 0.5;
    } else if(mean_x > mean_y) {
      pvalue = 0.0;
    } else {
      pvalue = 1.0;
    }
    return;
  }

  size_t v;
  datadefs::num_t sp,tvalue,ttrans;
    
  v = nreal_x + nreal_y - 2;
  sp = sqrt(((nreal_x-1) * var_x + (nreal_y-1) * var_y) / v);
  tvalue = (mean_x - mean_y) / (sp * sqrt(1.0 / nreal_x + 1.0 / nreal_y));
    
  ttrans = (tvalue+sqrt(pow(tvalue,2) + v)) / (2 * sqrt(pow(tvalue,2) + v));

  datadefs::regularized_betainc(ttrans,nreal_x - 1,pvalue);  
  pvalue = 1-pvalue;
    
}

//DEPRECATED

// !! Documentation: regularized beta function; not sure what the "inc" means
// !! in context, though. Is this part of an implementation of regularized beta
// !! functions, but lacking in additional functionality that exists outside of
// !! this method.

// !! Spurious deprecation: this isn't used in the existing codebase, but it
// !! has value within the API. Can we consider it "deprecated?"
void datadefs::regularized_betainc(const datadefs::num_t x,
                                   const size_t a,
                                   datadefs::num_t& ibval) {

  ibval = 0.0;
  
  datadefs::num_t jfac = 1;
  for(size_t i = 1; i < a; ++i) {
    jfac += log(static_cast<datadefs::num_t>(i));
  }
  
  datadefs::num_t kfac = 1;
  for(size_t i = a+1; i < 2*a; ++i) {
    kfac += log(static_cast<datadefs::num_t>(i));
  }
  
  for(size_t i = a; i < 2*a; ++i) {
    jfac += log(static_cast<datadefs::num_t>(i));
    kfac += log(static_cast<datadefs::num_t>(2*a - i));
    datadefs::num_t temp = kfac - jfac + i*log(x) + (2*a-1-i)*log(1-x);
    ibval += exp(temp);
      
    //cout << jfac << "\t" << kfac << "\t" << ibval << endl;
  }

  if(ibval > 1.0) {
    ibval = 1.0;
  }
  
}

// !! Documentation: performs a U-test across the input.

// !! Correctness: contains checks for "NAN" that alter the control flow of
// !!  this function. Consider a more meaningful sentinel and better handling
// !!  for NaN, as the presence of such greatly increases the returned p-value.

// !! Legibility: this function contains nested control flow that isn't immediately
// !! obvious to the reader. Consider refactoring.

void datadefs::utest(vector<datadefs::num_t> const& x,
                     vector<datadefs::num_t> const& y,
                     datadefs::num_t& pvalue) {
  
  num_t uvalue = 0.0;
  size_t m = 0;
  size_t n = 0;
  
  for(size_t i = 0; i < x.size(); ++i) {
    
    if (!datadefs::isNAN(x[i])) {
      ++m;
      for(size_t j = 0; j < y.size(); ++j) {
        
        if (!datadefs::isNAN(y[j])) {
          ++n;
          if (x[i] > y[j]) {
            uvalue += 1;
          }
        }
      }
    }
  }
  

  num_t mu = 1.0*m*n/2.0;
  num_t s = sqrt(1.0*m*n*(n+m+1)/12.0);

  //cout << uvalue << " " << mu << " " << s << endl;
  //cout << datadefs::erf( (uvalue-mu) / (s*sqrt(2.0)) )  << endl;
  
  pvalue = 1.0 - 0.5 * ( 1 + datadefs::erf( (uvalue-mu) / (s*sqrt(2.0)) ) );
  //cout << " ==> " << pvalue << endl;

}

// !! Documentation: computes the error function for the given value
datadefs::num_t datadefs::erf(datadefs::num_t x) {  

  // !! TODO: determine if raising an AssertionError or explicit exception is
  //     appropriate here.
  //assert(!datadefs::isNAN(x));
  //if (datadefs::isNAN(x)) {
  //  throw ERRNO_INVALID_ARGUMENT;
  //}
  
  num_t x2 = x*x;

  num_t sgn;
  if(x < 0.0) {
    sgn = -1.0;
  } else {
    sgn = 1.0;
  }

  return( sgn*sqrt(1.0 - exp(-x2*(4.0/datadefs::PI+datadefs::A*x2) / (1+datadefs::A*x2))) ); 

}

// !! Documentation: Spearman's rank correlation coefficient
// !! (http://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient)

// !! Correctness: this method isn't currently implemented. Consider a more
// !! meaningful exception than an assertion error, which will be silently
// !! compiled out if debugging is not enabled.
void datadefs::spearman_correlation(vector<datadefs::num_t> const& x,
                                    vector<datadefs::num_t> const& y,
                                    datadefs::num_t& corr) {
  assert(false);
}

// !! Documentation: Pearson product-moment correlation coefficient
// !! (http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient)

// !! Consistency: outputs NaN if x.size() == y.size() && y.size() == 0

void datadefs::pearson_correlation(vector<datadefs::num_t> const& x,
                                   vector<datadefs::num_t> const& y,
                                   datadefs::num_t& corr) {

  corr = 0.0;

  datadefs::num_t mu_x,se_x,mu_y,se_y;
  size_t nreal_x,nreal_y;

  vector<datadefs::num_t> x_real;
  vector<datadefs::num_t> y_real;

  size_t n = x.size();
  assert(n == y.size());
  
  for(size_t i = 0; i < n; ++i) {
    if(!datadefs::isNAN(x[i]) && !datadefs::isNAN(y[i])) {
      x_real.push_back(x[i]);
      y_real.push_back(y[i]);
    }
  }

  datadefs::sqerr(x_real,mu_x,se_x,nreal_x);
  datadefs::sqerr(y_real,mu_y,se_y,nreal_y);
  assert(nreal_x == nreal_y);
  
  for(size_t i = 0; i < nreal_x; ++i) {
    corr += ( x_real[i] - mu_x ) * ( y_real[i] - mu_y ); 
  }

  corr /= sqrt(se_x*se_y);

}



//DEPRECATED

// !! Documentation: computes the percentile for the value x, given an alpha
// !! score and "prc." I have no idea what "prc" means, here.

// !! Spurious deprecation: this isn't used in the existing codebase, but it
// !! has value within the API. Can we consider it "deprecated?"
void datadefs::percentile(vector<datadefs::num_t> x, 
                          const datadefs::num_t alpha, 
                          datadefs::num_t& prc) {

  sort(x.begin(),x.end());
  
  datadefs::num_t k((x.size()-1) * alpha);
  num_t f = floor(k);
  num_t c = ceil(k);
    
  if(fabs(f - c) < datadefs::EPS) {
    prc = x[static_cast<size_t>(k)];
  } else {
    num_t d0 = x[static_cast<size_t>(f)] * (c - k);
    num_t d1 = x[static_cast<size_t>(c)] * (k - f);
    prc = d0+d1;
  }
}
