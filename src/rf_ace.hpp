#ifndef RF_ACE_HPP
#define RF_ACE_HPP

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <set>
#include <string>
#include <vector>

#include "treedata.hpp"
#include "options.hpp"
#include "utils.hpp"

using namespace std;

namespace rface {

  void pruneFeatureSpace(Treedata& treeData, const options::General_options& gen_op) {

    size_t targetIdx = treeData.getFeatureIdx(gen_op.targetStr);

    if ( treeData.nRealSamples(targetIdx) == 0 ) {
      cerr << "Target feature '" << gen_op.targetStr << "' does not have any real samples!" << endl;
      exit(1);
    }
    
    if ( gen_op.whiteList != "" ) {
      
      cout << "===> Reading whitelist '" << gen_op.whiteList << "', please wait... " << flush;
      set<string> whiteFeatureNames = utils::readFeatureMask(treeData,gen_op.whiteList);
      cout << "DONE" << endl;
      cout << "===> Applying feature mask, removing " << treeData.nFeatures() - whiteFeatureNames.size()
	   << " / " << treeData.nFeatures() << " features, please wait... " << flush;
      treeData.whiteList(whiteFeatureNames);
      cout << "DONE" << endl;
    } 

    if ( gen_op.blackList != "" ) {
      
      cout << "===> Reading blacklist '" << gen_op.blackList << "', please wait... " << flush;
      set<string> blackFeatureNames = utils::readFeatureMask(treeData,gen_op.blackList);
      cout << "DONE" << endl;
      cout << "===> Applying blacklist, keeping " << treeData.nFeatures() - blackFeatureNames.size()
	   << " / " << treeData.nFeatures() << " features, please wait... " << flush;
      treeData.blackList(blackFeatureNames);
      cout << "DONE" << endl;
    }
    
    if ( gen_op.pruneFeatures ) {
      
      cout << "===> Pruning features with less than " << gen_op.pruneFeatures << " real samples... " << flush;
      size_t nFeaturesOld = treeData.nFeatures();
      utils::pruneFeatures(treeData,gen_op.targetStr,gen_op.pruneFeatures);
      cout << "DONE, " << nFeaturesOld - treeData.nFeatures() << " features ( "
	   << ( 100.0*(nFeaturesOld - treeData.nFeatures()) / nFeaturesOld ) << "% ) pruned" << endl;
      
    }

    if ( treeData.nFeatures() == 0 ) {
      cout << "All features were removed!" << endl;

      ofstream toLogFile(gen_op.log.c_str());

      toLogFile << "All features were removed!" << endl;

      toLogFile.close();

      exit(0);

    }
    
  }

  /*
    void validateRequiredParameters(const options::General_options& gen_op) {
    
    // Print help and exit if input file is not specified
    if ( gen_op.input == "" ) {
    cerr << "Input file not specified" << endl;
    options::printHelpHint();
    exit(1);
    }
    
    // Print help and exit if target index is not specified
    if ( !gen_op.isSet(gen_op.recombinePerms_s,gen_op.recombinePerms_l) && gen_op.targetStr == "" ) {
    cerr << "target not specified" << endl;
    options::printHelpHint();
    exit(1);
    }
    
    if ( gen_op.output == "" ) {
    cerr << "You forgot to specify an output file!" << endl;
    options::printHelpHint();
    exit(1);
    }
    
    }
  */

  void updateTargetStr(Treedata& treeData, options::General_options& gen_op) {

    // Check if the target is specified as an index
    int integer;
    if ( datadefs::isInteger(gen_op.targetStr,integer) ) {
      
      if ( integer < 0 || integer >= static_cast<int>( treeData.nFeatures() ) ) {
	cerr << "Feature index (" << integer << ") must be within bounds 0 ... " << treeData.nFeatures() - 1 << endl;
	exit(1);
      }
      
      // Extract the name of the feature, as upon masking the indices will become rearranged
      gen_op.targetStr = treeData.getFeatureName(static_cast<size_t>(integer));
      
    }

  }

  void updateMTry(Treedata& treeData, options::General_options& gen_op) {
    if ( gen_op.mTry == options::RF_DEFAULT_M_TRY ) {
      gen_op.mTry = static_cast<size_t>( floor(0.1*treeData.nFeatures()) );
    }
  }

  /*
    void printGeneralSetup(Treedata& treeData, const options::General_options& gen_op) {
    
    // After masking, it's safe to refer to features as indices
    // TODO: rf_ace.cpp: this should be made obsolete; instead of indices, use the feature headers
    size_t targetIdx = treeData.getFeatureIdx(gen_op.targetStr);
    
    size_t nAllFeatures = treeData.nFeatures();
    size_t nRealSamples = treeData.nRealSamples(targetIdx);
    num_t realFraction = 1.0*nRealSamples / treeData.nSamples();
    
    //Before number crunching, print values of parameters of RF-ACE
    cout << "General configuration:" << endl;
    cout << "    nfeatures" << setw(options::maxWidth-9) << "" << "= " << nAllFeatures << endl;
    cout << "    nsamples"  << setw(options::maxWidth-8) << "" << "= " << treeData.nRealSamples(targetIdx) << " / " << treeData.nSamples() << " ( " << 100.0 * ( 1 - realFraction ) << " % missing )" << endl;
    cout << "    tree type" << setw(options::maxWidth-9) << "" << "= ";
    if(treeData.isFeatureNumerical(targetIdx)) { cout << "Regression CART" << endl; } else { cout << treeData.nCategories(targetIdx) << "-class CART" << endl; }
    cout << "  --" << gen_op.dataDelimiter_l << setw( options::maxWidth - gen_op.dataDelimiter_l.size() ) << ""
    << "= '" << gen_op.dataDelimiter << "'" << endl;
    cout << "  --" << gen_op.headerDelimiter_l << setw( options::maxWidth - gen_op.headerDelimiter_l.size() ) << ""
    << "= '" << gen_op.headerDelimiter << "'" << endl;
    cout << "  --" << gen_op.input_l << setw( options::maxWidth - gen_op.input_l.size() ) << ""
    << "= " << gen_op.input << endl;
    cout << "  --" << gen_op.targetStr_l << setw( options::maxWidth - gen_op.targetStr_l.size() ) << ""
    << "= " << gen_op.targetStr << " ( index " << targetIdx << " )" << endl;
    cout << "  --" << gen_op.output_l << setw( options::maxWidth - gen_op.output_l.size() ) << ""
    << "= "; if ( gen_op.output != "" ) { cout << gen_op.output << endl; } else { cout << "NOT SET" << endl; }
    cout << "  --" << gen_op.log_l << setw( options::maxWidth - gen_op.log_l.size() ) << ""
    << "= "; if( gen_op.log != "" ) { cout << gen_op.log << endl; } else { cout << "NOT SET" << endl; }
    cout << "  --" << gen_op.seed_l << setw( options::maxWidth - gen_op.seed_l.size() ) << ""
    << "= " << gen_op.seed << endl;
    cout << endl;
    
    cout << "Stochastic Forest configuration:" << endl;
    cout << "  --" << gen_op.nTrees_l << setw( options::maxWidth - gen_op.nTrees_l.size() ) << ""
    << "= "; if(gen_op.nTrees == 0) { cout << "DEFAULT" << endl; } else { cout << gen_op.nTrees << endl; }
    cout << "  --" << gen_op.mTry_l << setw( options::maxWidth - gen_op.mTry_l.size() ) << ""
    << "= " << gen_op.mTry << endl;
    cout << "  --" << gen_op.nMaxLeaves_l << setw( options::maxWidth - gen_op.nMaxLeaves_l.size() ) << ""
    << "= " << gen_op.nMaxLeaves << endl;
    cout << "  --" << gen_op.nodeSize_l << setw( options::maxWidth - gen_op.nodeSize_l.size() ) << ""
    << "= "; if(gen_op.nodeSize == 0) { cout << "DEFAULT" << endl; } else { cout << gen_op.nodeSize << endl; }
    cout << "  --" << gen_op.shrinkage_l << setw( options::maxWidth - gen_op.shrinkage_l.size() ) << ""
    << "= " << gen_op.shrinkage << endl;
    cout << endl;
    
    }
  */
  
}

#endif
