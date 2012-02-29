#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <set>
#include <string>
#include <istream>
//#include <algorithm>
#include <cstdlib>
#include "datadefs.hpp"
#include "treedata.hpp"

using namespace std;
using datadefs::num_t;


namespace utils {

  // Removes missing values from the provided data vector
  vector<num_t> removeNANs(vector<num_t> x);

  int str2int(const string& str);

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

  string join(const vector<string>& items, const char delimiter);

  // WILL BECOME OBSOLETE
  set<string> readFeatureMask(Treedata& treeData, const string& fileName);

  // MAY BECOME OBSOLETE
  void pruneFeatures(Treedata& treeData, const string& targetName, const size_t minSamples);

  void filterSort(const bool isIncreasingOrder,
		  vector<num_t>& data,
		  vector<size_t>& refIcs);
  
}

#endif
