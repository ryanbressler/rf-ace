#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <set>
#include <map>
#include <string>
#include <istream>
#include <sstream>
#include <cstdlib>
#include <unordered_set>
#include "datadefs.hpp"
//#include "math.hpp"
#include "distributions.hpp"

using namespace std;
using datadefs::num_t;

class Treedata;

namespace utils {

  string tolower(const string& str);

  string suffix(const string& str);

  // Removes missing values from the provided data vector
  vector<num_t> removeNANs(vector<num_t> x);
    
  // Chomps a string, i.e. removes all the trailing end-of-line characters
  string chomp(const string& str, const string& eof = "\r\n");

  // Removes leading and trailing whitespace characters
  string trim(const string& str, const string& wh = " ");

  unordered_set<string> keys(const string& str, const char delimiter);
  
  // A sophisticated parser that extracts a key-value pairs from a string
  map<string,string> parse(const string& str, 
			   const char delimiter, 
			   const char separator, 
			   const char comment);

  map<string,string> parse(istream& streamObj,
			   const char delimiter,
			   const char separator,
			   const char comment);

  unordered_set<uint32_t> hashText(const string& text);
  
  // Splits a delimited string
  vector<string> split(const string& str, const char delimiter, const string& wh = " ");

  // Splits a delimited stream
  vector<string> split(istream& streamObj, const char delimiter, const string& wh = " ");

  // Reads a list of items from a file
  vector<string> readListFromFile(const string& fileName, const char delimiter);

  template<typename StartIterator, typename StopIterator>
  inline void write(ostream& os, StartIterator startIt, StopIterator stopIt, const char delimiter = ' ') {

    if ( startIt != stopIt ) {
      os << *startIt;
      ++startIt;
    }
    
    while ( startIt != stopIt ) {
      os << delimiter << *startIt;
      ++startIt;
    }
  }

  void filterSort(const bool isIncreasingOrder,
		  vector<num_t>& data,
		  vector<size_t>& refIcs);
  
  string num2str(const num_t x);

  void strv2numv(const vector<string>& strvec,
		 vector<datadefs::num_t>& numvec);

  void strv2catv(const vector<string>& strvec, 
		 vector<datadefs::num_t>& catvec, 
		 map<string,datadefs::num_t>& mapping, 
		 map<datadefs::num_t,string>& backMapping);

  void sortDataAndMakeRef(const bool isIncreasingOrder,
			  vector<datadefs::num_t>& data,
			  vector<size_t>& refIcs);

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
    
  template <typename T>
  T str2(const string& str) {

    if( datadefs::isNAN_STR(str) ) {
      return( static_cast<T>(datadefs::NUM_NAN) );
    }
    
    stringstream ss( chomp(str) );
    T ret;
    ss >> ret;
    
    if ( ss.fail() || ss.bad() || !ss.eof() ) {
      cerr << "utils::convert::str2<T>() -- input '" << str
	   << "' incorrectly formatted for conversion to type T" << endl;
      exit(1);
    }
    
    return( ret );
  }

  vector<size_t> range(const size_t n);

  istream& safeGetline(istream& is, string& t);

  vector<vector<size_t> > splitRange(const size_t nElements, const size_t nSplits);

  template<typename T>
  void permute(vector<T>& data, distributions::Random* random) {
    
    // Permute indices
    for (size_t i = 0; i < data.size(); ++i) {
      size_t j = random->integer() % (i + 1);
      T temp = data[i];
      data[i] = data[j];
      data[j] = temp;
    }

  }

  num_t numericalFeatureSplitsNumericalTarget(const vector<num_t>& tv,
					      const vector<num_t>& fv,
					      const size_t minSamples,
					      size_t& splitIdx);
  
  num_t numericalFeatureSplitsCategoricalTarget(const vector<num_t>& tv,
						const vector<num_t>& fv,
						const size_t minSamples,
						size_t& splitIdx);
  
  num_t categoricalFeatureSplitsNumericalTarget(const vector<num_t>& tv,
						const vector<num_t>& fv,
						const size_t minSamples,
						map<num_t,vector<size_t> >& fmap_left,
						map<num_t,vector<size_t> >& fmap_right);
 
  num_t categoricalFeatureSplitsNumericalTarget2(const vector<num_t>& tv,
						 const vector<num_t>& fv,
						 const size_t minSamples,
						 const vector<num_t>& catOrder,
						 map<num_t,vector<size_t> >& fmap_left,
						 map<num_t,vector<size_t> >& fmap_right);
 
  num_t categoricalFeatureSplitsCategoricalTarget(const vector<num_t>& tv,
						  const vector<num_t>& fv,
						  const size_t minSamples,
						  map<num_t,vector<size_t> >& fmap_left,
						  map<num_t,vector<size_t> >& fmap_right);

  num_t categoricalFeatureSplitsCategoricalTarget2(const vector<num_t>& tv,
						   const vector<num_t>& fv,
						   const size_t minSamples,
						   const vector<num_t>& catOrder,
						   map<num_t,vector<size_t> >& fmap_left,
						   map<num_t,vector<size_t> >& fmap_right);
  
  
}

#endif
