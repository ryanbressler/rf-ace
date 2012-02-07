#include <sstream>

#include "utils.hpp"

vector<num_t> utils::removeNANs(const vector<num_t>& x) {

  vector<num_t> trimmed(x.size());

  size_t nRetained = 0;
  for(size_t i = 0; i < x.size(); ++i) {
    if( !datadefs::isNAN(x[i]) ) {
      trimmed[nRetained] = x[i];
      ++nRetained;
    }
  }
  trimmed.resize(nRetained);
  if ( nRetained == 0 ) {
    cout << "utils::removeNANs() -- data has 0 real values!" << endl;
  }

  return(trimmed);
}

// Makes a copy of the string, chomps, and returns
string utils::chomp(const string& str) {
  
  // Chop at the first newline character, if it exists
  int crIdx = str.find("\r");
  int lfIdx = str.find("\n");
  int terminatorIdx = crIdx;
  if (lfIdx != -1 && lfIdx < crIdx) {
    terminatorIdx = lfIdx;
  }
  
  string ret(str);
  
  return( terminatorIdx != -1 ? ret.substr(0,terminatorIdx) : ret );
}

vector<string> utils::readListFromFile(const string& fileName, const char delimiter) {
  
  ifstream listStream;
  listStream.open(fileName.c_str());
  assert(listStream.good());
  
  string newItem("");
  vector<string> items;
  
  while ( getline(listStream,newItem,delimiter) ) {
    
    newItem = utils::chomp(newItem);
    
    items.push_back(newItem);
    
  }

  return( items );
  
}

set<string> utils::readFeatureMask(Treedata& treeData, const string& fileName) {
  
  ifstream featurestream;
  featurestream.open(fileName.c_str());
  assert(featurestream.good());
  
  string newFeature;
  
  set<string> featureMaskSet;
  
  getline(featurestream,newFeature);

  int integer;
  if ( datadefs::isInteger(newFeature,integer) ) {

    size_t foo = static_cast<size_t>(integer);
    newFeature = treeData.getFeatureName(foo);

    featureMaskSet.insert( newFeature );

    while ( getline(featurestream,newFeature) ) {

      stringstream ss(newFeature);
      ss >> foo;
      featureMaskSet.insert( treeData.getFeatureName(foo) );
    }

  } else {

    string newFeature;

    featureMaskSet.insert( newFeature );

    while ( getline(featurestream,newFeature) ) {
      featureMaskSet.insert(newFeature);
    }
  }

  return( featureMaskSet );
}

void utils::pruneFeatures(Treedata& treeData, const string& targetName, const size_t minSamples) {

  set<string> featureNames;

  size_t targetIdx = treeData.getFeatureIdx(targetName);

  for ( size_t featureIdx = 0; featureIdx < treeData.nFeatures(); ++featureIdx ) {

    if ( treeData.nRealSamples(targetIdx,featureIdx) < minSamples ) {
      featureNames.insert(treeData.getFeatureName(featureIdx));
    }

  }

  treeData.blackList(featureNames);

}



