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

int utils::str2int(const string& str) {

  stringstream ss( utils::chomp(str) );
  int ret;
  ss >> ret;

  return( ret );

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

set<string> utils::keys(const string& str, const char delimiter) {

  set<string> ret;

  vector<string> items = utils::split(str,delimiter);

  for ( size_t i = 0; i < items.size(); ++i ) {
    ret.insert(items[i]);
  }
  
  return( ret );

}

/*
  map<string,string> utils::keys2vals(const string& str,
  const char delimiter,
  const char separator) {
  
  map<string,string> ret;
  
  vector<string> items = utils::split(str,delimiter);
  
  for( size_t i = 0; i < items.size(); ++i ) {
  vector<string> foo = utils::split(items[i],separator);
  ret[foo[0]] = foo[1];
  }
  
  return( ret );
  
  }
*/

map<string,string> utils::parse(const string& str,
				const char delimiter,
				const char separator,
				const char comment) {

  stringstream streamObj(str);

  return( utils::parse(streamObj,delimiter,separator,comment) );

}

map<string,string> utils::parse(istream& streamObj,
				const char delimiter,
				const char separator,
				const char comment) {

  map<string,string> ret;

  string key;

  while ( !streamObj.eof() ) {

    // Parse the key
    getline(streamObj,key,separator);
    assert( streamObj.good() );
    assert( !streamObj.eof() );
    
    // Peek the next characeter and check if it's a comment 
    if ( streamObj.peek() == comment ) {

      // ignore the comment character...
      streamObj.ignore();

      // ... and get the value for the key
      getline(streamObj,ret[key],comment);

      assert( ret.find(key) != ret.end() );

      // If the next character is a delimiter, ignore it
      if ( streamObj.peek() == delimiter ) {
	streamObj.ignore();
      } 

    } else {
      
      // Without the comment character we just read until the next delimiter
      getline(streamObj,ret[key],delimiter);
      
    }
 
  }
  
  // The possible carriage return and end-of-line characters need to be removed
  ret[key] = utils::chomp(ret[key]);
    
  return( ret );

}

vector<string> utils::split(const string& str, const char delimiter) {
  stringstream streamObj(str);
  return( utils::split(streamObj,delimiter) );
}

vector<string> utils::split(istream& streamObj, const char delimiter) {

  string newItem("");
  vector<string> items;

  while ( getline(streamObj,newItem,delimiter) ) {
    newItem = utils::chomp(newItem);
    items.push_back(newItem);
  }

  return( items );

}

vector<string> utils::readListFromFile(const string& fileName, const char delimiter) {
  
  ifstream streamObj( fileName.c_str() );
  assert(streamObj.good());
  
  return( utils::split(streamObj,delimiter) );
}


string utils::join(const vector<string>& items, const char delimiter) {

  string ret("");

  if ( items.size() == 0 ) {
    return( ret );
  }

  ret = items[0];

  for( size_t i = 1; i < items.size(); ++i ) {
    ret.append(delimiter+items[i]);
  }

  return( ret );

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

    //string newFeature;

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

void utils::filterSort(const bool isIncreasingOrder,
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
  
  pairedv.erase(remove_if(pairedv.begin(),pairedv.end(),&datadefs::pairedIsNAN), pairedv.end());
  
  if(isIncreasingOrder) {
    sort(pairedv.begin(),pairedv.end(),datadefs::increasingOrder<size_t>());
  } else {
    sort(pairedv.begin(),pairedv.end(),datadefs::decreasingOrder<size_t>());
  }
  
  datadefs::separate_pairedv<num_t,size_t>(pairedv,data,refIcs);
}



