#include <vector>
#include <set>
#include <string>
#include "datadefs.hpp"
#include "treedata.hpp"

using namespace std;
using datadefs::num_t;


namespace utils {

  // Removes missing values from the provided data vector
  vector<num_t> removeNANs(const vector<num_t>& x);

  // Chomps a string, i.e. removes all the trailing end-of-line characters
  string chomp(const string& str);

  // Reads a delimited list of strings from a file
  vector<string> readListFromFile(const string& fileName, const char delimiter);

  set<string> readFeatureMask(Treedata& treeData, const string& fileName);

  void pruneFeatures(Treedata& treeData, const string& targetName, const size_t minSamples);

  

}

