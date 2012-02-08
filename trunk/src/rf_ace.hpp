#ifndef RF_ACE_HPP
#define RF_ACE_HPP

#include <cstdlib>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "treedata.hpp"
#include "options.hpp"
#include "utils.hpp"

using namespace std;

namespace rface {

  void pruneFeatureSpace(Treedata& treeData, const options::General_options& gen_op) {
    
    if ( gen_op.whiteList != "" ) {
      
      cout << "Reading whitelist '" << gen_op.whiteList << "', please wait... " << flush;
      set<string> whiteFeatureNames = utils::readFeatureMask(treeData,gen_op.whiteList);
      cout << "DONE" << endl;
      cout << "Applying feature mask, removing " << treeData.nFeatures() - whiteFeatureNames.size()
	   << " / " << treeData.nFeatures() << " features, please wait... " << flush;
      treeData.whiteList(whiteFeatureNames);
      cout << "DONE" << endl;
    } 
    
    if ( gen_op.blackList != "" ) {
      
      cout << "Reading blacklist '" << gen_op.blackList << "', please wait... " << flush;
      set<string> blackFeatureNames = utils::readFeatureMask(treeData,gen_op.blackList);
      cout << "DONE" << endl;
      cout << "Applying blacklist, keeping " << treeData.nFeatures() - blackFeatureNames.size()
	   << " / " << treeData.nFeatures() << " features, please wait... " << flush;
      treeData.blackList(blackFeatureNames);
      cout << "DONE" << endl;
    }
    
    if ( gen_op.pruneFeatures ) {
      
      cout << "Pruning features with less than " << gen_op.pruneFeatures << " real samples... " << flush;
      size_t nFeaturesOld = treeData.nFeatures();
      utils::pruneFeatures(treeData,gen_op.targetStr,gen_op.pruneFeatures);
      cout << "DONE, " << nFeaturesOld - treeData.nFeatures() << " features ( "
	   << ( 100.0*(nFeaturesOld - treeData.nFeatures()) / nFeaturesOld ) << "% ) pruned" << endl;
      
    }
    
  }
  
  
}

#endif
