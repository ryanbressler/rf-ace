#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <set>
#include <string>
#include <istream>
#include <sstream>
#include <cstdlib>
#include "datadefs.hpp"
#include "treedata.hpp"

using namespace std;
using datadefs::num_t;


namespace utils {

  // Generate seeds for random number generators
  int generateSeed();
  
  // Removes missing values from the provided data vector
  vector<num_t> removeNANs(vector<num_t> x);

  
    
  // Chomps a string, i.e. removes all the trailing end-of-line characters
  string chomp(const string& str);

  set<string> keys(const string& str, const char delimiter);
  
  // A sophisticated parser that extracts a key-value pairs from a string
  map<string,string> parse(const string& str, 
			   const char delimiter, 
			   const char separator, 
			   const char comment);

  map<string,string> parse(istream& streamObj,
			   const char delimiter,
			   const char separator,
			   const char comment);
  
  // Splits a delimited string
  vector<string> split(const string& str, const char delimiter);

  // Splits a delimited stream
  vector<string> split(istream& streamObj, const char delimiter);

  // Reads a list of items from a file
  vector<string> readListFromFile(const string& fileName, const char delimiter);

  template <typename InputIterator>
  string join(InputIterator begin, InputIterator end, const char delimiter) {

    string ret("");
    
    if ( begin == end ) {
      return(ret);
    }

    stringstream ss;
    ss << *begin;
    ++begin;
    
    while( begin != end ) {
      
      ss << delimiter << *begin;

      ++begin;

    }

    ss >> ret;
    return( ret );

  }

  // WILL BECOME OBSOLETE
  set<string> readFeatureMask(Treedata& treeData, const string& fileName);

  // MAY BECOME OBSOLETE
  void pruneFeatures(Treedata& treeData, const string& targetName, const size_t minSamples);

  void filterSort(const bool isIncreasingOrder,
		  vector<num_t>& data,
		  vector<size_t>& refIcs);
  
  string num2str(const num_t x);

  template <typename T>
  T str2(const string& str) {
    
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
  
}

#endif
