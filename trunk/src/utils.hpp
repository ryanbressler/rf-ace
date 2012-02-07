#include <vector>
#include <set>
#include <string>
#include "datadefs.hpp"
#include "treedata.hpp"

using namespace std;
using datadefs::num_t;


namespace utils {

  vector<num_t> removeNANs(const vector<num_t>& x);

  set<string> readFeatureMask(Treedata& treeData, const string& fileName);

  void pruneFeatures(Treedata& treeData, const string& targetName, const size_t minSamples);

}

