#ifndef RF_ACE_HPP
#define RF_ACE_HPP

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <set>
#include <string>
#include <vector>
#include <stdio.h>
#include <ctime>

#include "stochasticforest.hpp"
#include "treedata.hpp"
#include "options.hpp"
#include "utils.hpp"
#include "math.hpp"

using namespace std;

namespace rface {

  void printHeader(ostream& out) {
    out << endl
	<< "-----------------------------------------------------------" << endl
	<< "|  RF-ACE version:  1.0.7, Aug 28 2012                    |" << endl
	<< "|    Compile date:  " << __DATE__ << ", " << __TIME__ << "                 |" << endl
	<< "|   Report issues:  code.google.com/p/rf-ace/issues/list  |" << endl
	<< "-----------------------------------------------------------" << endl
	<< endl;
  }

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
      
      // Add the target feature into the white list, otherwise it may get removed
      whiteFeatureNames.insert(gen_op.targetStr);

      treeData.whiteList(whiteFeatureNames);
      cout << "DONE" << endl;
    } 

    if ( gen_op.blackList != "" ) {
      
      cout << "===> Reading blacklist '" << gen_op.blackList << "', please wait... " << flush;
      set<string> blackFeatureNames = utils::readFeatureMask(treeData,gen_op.blackList);
      cout << "DONE" << endl;
      cout << "===> Applying blacklist, keeping " << treeData.nFeatures() - blackFeatureNames.size()
	   << " / " << treeData.nFeatures() << " features, please wait... " << flush;

      // Remove the target feature from the black list, otherwise it will get removed
      if ( blackFeatureNames.find(gen_op.targetStr) != blackFeatureNames.end() ) {
	cout << " Target found in the blacklist -- omitting... " << flush;
	blackFeatureNames.erase(gen_op.targetStr);
      }
      
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

  void printGeneralSetup(Treedata& treeData, const options::General_options& gen_op) {

    // After masking, it's safe to refer to features as indices
    // TODO: rf_ace.cpp: this should be made obsolete; instead of indices, use the feature headers
    size_t targetIdx = treeData.getFeatureIdx(gen_op.targetStr);

    size_t nAllFeatures = treeData.nFeatures();
    size_t nRealSamples = treeData.nRealSamples(targetIdx);
    num_t realFraction = 1.0*nRealSamples / treeData.nSamples();

    //Before number crunching, print values of parameters of RF-ACE
    cout << endl;
    cout << "Input data:" << endl;
    cout << " - " << nAllFeatures << " features" << endl;
    cout << " - " << treeData.nRealSamples(targetIdx) << " samples / " << treeData.nSamples() << " ( " << 100.0 * ( 1 - realFraction ) << " % missing )" << endl;

  }

  void setEnforcedForestParameters(Treedata& treeData, options::General_options& gen_op) {

    if ( gen_op.modelType == options::RF ) {

      size_t targetIdx = treeData.getFeatureIdx(gen_op.targetStr);

      // Allow trees to grow to maximal depth, if not told otherwise
      gen_op.setIfNotSet(gen_op.nMaxLeaves_s,gen_op.nMaxLeaves_l,gen_op.nMaxLeaves,treeData.nRealSamples(targetIdx));

      // RF mTry is by default set to 10% of features
      gen_op.setIfNotSet(gen_op.mTry_s,gen_op.mTry_l,gen_op.mTry,static_cast<size_t>(0.1*treeData.nFeatures()));

      // Minimum mTry is 1
      if ( gen_op.mTry < 1 ) {
	gen_op.mTry = 1;
      }

    } else if ( gen_op.modelType == options::CART ) {

      // In CART mode only one tree is grown
      gen_op.nTrees = 1;
    }

  }

  StochasticForest buildPredictor(Treedata& trainData, options::General_options& gen_op) {
    
    updateTargetStr(trainData,gen_op);
    
    pruneFeatureSpace(trainData,gen_op);
    
    setEnforcedForestParameters(trainData,gen_op);
    
    // We never want to use contrasts when we are building a predictor
    //gen_op.useContrasts = false;
    
    printGeneralSetup(trainData,gen_op);
    
    gen_op.print();
    
    gen_op.validateParameters();
    
    if ( gen_op.modelType == options::RF ) {
      cout << "===> Growing RF predictor... " << flush;
    } else if ( gen_op.modelType == options::GBT ) {
      cout << "===> Growing GBT predictor... " << flush;
    } else if ( gen_op.modelType == options::CART ) {
      cout << "===> Growing CART predictor... " << flush;
    } else {
      cerr << "Unknown forest type!" << endl;
      exit(1);
    }
    
    StochasticForest SF(&trainData,gen_op);
    cout << "DONE" << endl << endl;
    
    if ( gen_op.modelType == options::GBT ) {
      cout << "GBT diagnostics disabled temporarily" << endl << endl;
      return( SF );
    }
    
    size_t targetIdx = trainData.getFeatureIdx(gen_op.targetStr);
    vector<num_t> data = utils::removeNANs(trainData.getFeatureData(targetIdx));

    num_t oobError = SF.getOobError();
    num_t ibOobError =  SF.getError();
    
    cout << "RF training error measures (NULL == no model):" << endl;
    if ( trainData.isFeatureNumerical(targetIdx) ) {
      num_t nullError = math::var(data);
      cout << "              NULL std = " << sqrt( nullError ) << endl;
      cout << "               OOB std = " << sqrt( oobError ) << endl;
      cout << "            IB+OOB std = " << sqrt( ibOobError ) << endl;
      cout << "  % explained by model = " << 1 - oobError / nullError << " = 1 - (OOB var) / (NULL var)" << endl;
    } else {
      num_t nullError = math::nMismatches( data, math::mode(data) );
      cout << "       NULL % mispred. = " << 1.0 * nullError / data.size() << endl;
      cout << "        OOB % mispred. = " << oobError << endl;
      cout << "     IB+OOB % mispred. = " << ibOobError << endl;
      cout << "  % explained by model = " << 1 - oobError / nullError << " ( 1 - (OOB # mispred.) / (NULL # mispred.) )" << endl;
    }
    cout << endl;

    return( SF );

  }


  
}

#endif
