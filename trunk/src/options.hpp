#ifndef OPTIONS_HPP
#define OPTIONS_HPP

#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include "argparse.hpp"
#include "datadefs.hpp"
#include "treedata.hpp"
#include "utils.hpp"

using namespace std;
using datadefs::num_t;

namespace options {

  void printHeader(ostream& output) {

    output << endl;
    output << "-----------------------------------------------------" << endl;
    output << "|  RF-ACE version:  1.0.5, April 24th, 2012         |" << endl;
    output << "|    Project page:  http://code.google.com/p/rf-ace |" << endl;
    output << "|     Report bugs:  timo.p.erkkila@tut.fi           |" << endl;
    output << "-----------------------------------------------------" << endl;
    output << endl;
  }

  void printHelpHint() {
    cout << endl;
    cout << "To get started, type \"-h\" or \"--help\"" << endl;
  }
  
  // Default general configuration
  const bool   GENERAL_DEFAULT_PRINT_HELP = false;
  const bool   GENERAL_DEFAULT_IS_FILTER = false;
  const char   GENERAL_DEFAULT_DATA_DELIMITER = '\t';
  const char   GENERAL_DEFAULT_HEADER_DELIMITER = ':';
  const size_t GENERAL_DEFAULT_MIN_SAMPLES = 5;
  const int    GENERAL_DEFAULT_SEED = -1;

  // Statistical test default configuration
  const size_t ST_DEFAULT_N_PERMS = 20;
  const num_t  ST_DEFAULT_P_VALUE_THRESHOLD = 0.05;
  const bool   ST_DEFAULT_NO_SORT = false;

  // Random Forest default configuration
  const size_t RF_DEFAULT_N_TREES = 100;
  const size_t RF_DEFAULT_M_TRY = 0;
  const size_t RF_DEFAULT_N_MAX_LEAVES = 100;
  const size_t RF_DEFAULT_NODE_SIZE = 3;
  const num_t  RF_DEFAULT_IN_BOX_FRACTION = 1.0;
  const num_t  RF_DEFAULT_SAMPLE_WITH_REPLACEMENT = true;
  const bool   RF_DEFAULT_USE_CONTRASTS = true;
  const bool   RF_DEFAULT_IS_RANDOM_SPLIT = true;
  const num_t  RF_DEFAULT_SHRINKAGE = 0.0;

  // Gradient Boosting Trees default configuration
  const size_t GBT_DEFAULT_N_TREES = 100;
  const size_t GBT_DEFAULT_M_TRY = 0;
  const size_t GBT_DEFAULT_N_MAX_LEAVES = 6;
  const size_t GBT_DEFAULT_NODE_SIZE = 3;
  const num_t  GBT_DEFAULT_IN_BOX_FRACTION = 0.5;
  const num_t  GBT_DEFAULT_SAMPLE_WITH_REPLACEMENT = false;
  const bool   GBT_DEFAULT_USE_CONTRASTS = false;
  const bool   GBT_DEFAULT_IS_RANDOM_SPLIT = false;
  const num_t  GBT_DEFAULT_SHRINKAGE = 0.1;

  // Determines the amount of indentation in help print-outs
  const size_t maxWidth = 20;
  
  struct General_options {

    // I/O related and general paramters
    bool printHelp; const string printHelp_s; const string printHelp_l;
    bool isFilter; const string isFilter_s; const string isFilter_l;
    string input; const string input_s; const string input_l;
    string output; const string output_s; const string output_l;
    string targetStr; const string targetStr_s; const string targetStr_l;
    string whiteList; const string whiteList_s; const string whiteList_l;
    string blackList; const string blackList_s; const string blackList_l;
    string predictionData; const string predictionData_s; const string predictionData_l;
    string log; const string log_s; const string log_l;
    char dataDelimiter; const string dataDelimiter_s; const string dataDelimiter_l;
    char headerDelimiter; const string headerDelimiter_s; const string headerDelimiter_l;
    size_t pruneFeatures; const string pruneFeatures_s; const string pruneFeatures_l;
    int seed; string seed_s; string seed_l;

    // Random Forest related parameters
    size_t  nTrees; const string nTrees_s; const string nTrees_l;
    size_t  mTry; const string mTry_s; const string mTry_l;
    size_t  nMaxLeaves; const string nMaxLeaves_s; const string nMaxLeaves_l;
    size_t  nodeSize; const string nodeSize_s; const string nodeSize_l;
    num_t   shrinkage; const string shrinkage_s; const string shrinkage_l;

    // Statistical test related parameters
    size_t nPerms; const string nPerms_s; const string nPerms_l;
    num_t pValueThreshold; const string pValueThreshold_s; const string pValueThreshold_l;
    bool noSort; const string noSort_s; const string noSort_l;
    string contrastOutput; const string contrastOutput_s; const string contrastOutput_l;

    General_options():
      printHelp(GENERAL_DEFAULT_PRINT_HELP),printHelp_s("h"),printHelp_l("help"),
      isFilter(GENERAL_DEFAULT_IS_FILTER),isFilter_s("L"),isFilter_l("filter"),
      input(""),input_s("I"),input_l("input"),
      output(""),output_s("O"),output_l("output"),
      targetStr(""),targetStr_s("i"),targetStr_l("target"),
      whiteList(""),whiteList_s("W"),whiteList_l("whitelist"),
      blackList(""),blackList_s("B"),blackList_l("blacklist"),
      predictionData(""),predictionData_s("T"),predictionData_l("test"),
      log(""),log_s("L"),log_l("log"),
      dataDelimiter(GENERAL_DEFAULT_DATA_DELIMITER),dataDelimiter_s("D"),dataDelimiter_l("data_delim"),
      headerDelimiter(GENERAL_DEFAULT_HEADER_DELIMITER),headerDelimiter_s("H"),headerDelimiter_l("head_delim"),
      pruneFeatures(GENERAL_DEFAULT_MIN_SAMPLES),pruneFeatures_s("X"),pruneFeatures_l("prune_features"),
      seed(GENERAL_DEFAULT_SEED),seed_s("S"),seed_l("seed"),
      // Random Forest related parameters
      nTrees_s("n"),nTrees_l("ntrees"),
      mTry_s("m"),mTry_l("mtry"),
      nMaxLeaves_s("a"),nMaxLeaves_l("nmaxleaves"),
      nodeSize_s("s"),nodeSize_l("nodesize"),
      shrinkage_s("k"),shrinkage_l("shrinkage"),
      // Statistical test related parameters
      nPerms(ST_DEFAULT_N_PERMS),nPerms_s("p"),nPerms_l("nperms"),
      pValueThreshold(ST_DEFAULT_P_VALUE_THRESHOLD),pValueThreshold_s("t"),pValueThreshold_l("pthreshold"), 
      noSort(ST_DEFAULT_NO_SORT),noSort_s("N"),noSort_l("noSort"),
      contrastOutput(""),contrastOutput_s("C"),contrastOutput_l("contrast_output") { setRFDefaults(); }
    
    void loadUserParams(ArgParse& parser) {
            
      // I/O related and general parameters
      parser.getFlag(printHelp_s, printHelp_l, printHelp);
      parser.getFlag(isFilter_s,isFilter_l,isFilter);
      parser.getArgument<string>(input_s, input_l, input);
      parser.getArgument<string>(targetStr_s, targetStr_l, targetStr);
      parser.getArgument<string>(output_s, output_l, output);
      parser.getArgument<string>(whiteList_s, whiteList_l, whiteList);
      parser.getArgument<string>(blackList_s, blackList_l, blackList);
      parser.getArgument<string>(predictionData_s, predictionData_l, predictionData);
      parser.getArgument<string>(log_s,log_l,log);
      string dataDelimiter,headerDelimiter;
      parser.getArgument<string>(dataDelimiter_s, dataDelimiter_l, dataDelimiter);
      parser.getArgument<string>(headerDelimiter_s, headerDelimiter_l, headerDelimiter);
      parser.getArgument<size_t>(pruneFeatures_s, pruneFeatures_l, pruneFeatures);
      stringstream ss(dataDelimiter);
      ss >> dataDelimiter;

      ss.clear();
      ss.str("");
      ss << headerDelimiter;
      ss >> headerDelimiter;

      parser.getArgument<int>(seed_s, seed_l, seed);
      
      // If no seed was provided, generate one
      if ( seed == GENERAL_DEFAULT_SEED ) {
	seed = utils::generateSeed();
      }
    
      // Random Forest related parameters
      parser.getArgument<size_t>( nTrees_s, nTrees_l, nTrees );
      parser.getArgument<size_t>( mTry_s, mTry_l, mTry );
      parser.getArgument<size_t>( nMaxLeaves_s, nMaxLeaves_l, nMaxLeaves );
      parser.getArgument<size_t>( nodeSize_s, nodeSize_l, nodeSize );
      parser.getArgument<num_t>(  shrinkage_s, shrinkage_l, shrinkage );

      // Statistical test parameters
      parser.getArgument<size_t>(nPerms_s, nPerms_l, nPerms);
      parser.getArgument<num_t>(pValueThreshold_s,  pValueThreshold_l,  pValueThreshold);
      parser.getFlag(noSort_s, noSort_l, noSort);
      parser.getArgument<string>(contrastOutput_s, contrastOutput_l, contrastOutput);

      if ( nPerms > 1 && nPerms < 6 ) {
        cerr << "Use more than 5 permutations in statistical test!" << endl;
        exit(1);
      }

      if ( pValueThreshold < 0.0 || pValueThreshold > 1.0 ) {
        cerr << "P-value threshold in statistical test must be within 0...1" << endl;
        exit(1);
      }

  
    }

    void help() {

      cout << "GENERAL ARGUMENTS:" << endl;
      cout << " -" << input_s << " / --" << input_l << setw( maxWidth - input_l.size() )
           << " " << "Input data file (.afm or .arff) or predictor forest file (.sf)" << endl;
      cout << " -" << targetStr_s << " / --" << targetStr_l << setw( maxWidth - targetStr_l.size() )
           << " " << "Target, specified as integer or string that is to be matched with the content of input" << endl;
      cout << " -" << output_s << " / --" << output_l << setw( maxWidth - output_l.size() )
           << " " << "Output file" << endl;
      cout << " -" << predictionData_s << " / --" << predictionData_l << setw( maxWidth - predictionData_l.size() )
	   << " " << "[Prediction only] Test data file (.afm or .arff) for prediction" << endl;
      cout << " -" << log_s << " / --" << log_l << setw( maxWidth - log_l.size() )
           << " " << "Log output file" << endl;
      cout << " -" << whiteList_s << " / --" << whiteList_l << setw( maxWidth - whiteList_l.size() )
           << " " << "White list of features to be included in the input file(s)." << endl;
      cout << " -" << blackList_s << " / --" << blackList_l << setw( maxWidth - blackList_l.size() )
           << " " << "Black list of features to be excluded from the input file(s)." << endl;
      cout << " -" << dataDelimiter_s << " / --" << dataDelimiter_l << setw( maxWidth - dataDelimiter_l.size() )
           << " " << "[AFM only] Data delimiter (default \\t)" << endl;
      cout << " -" << headerDelimiter_s << " / --" << headerDelimiter_l << setw( maxWidth - headerDelimiter_l.size() )
           << " " << "[AFM only] Header delimiter that separates the N and C symbols in feature headers from the rest (default " << GENERAL_DEFAULT_HEADER_DELIMITER << ")" << endl;
      cout << " -" << pruneFeatures_s << " / --" << pruneFeatures_l << setw( maxWidth - pruneFeatures_l.size() )
           << " " << "Features with less than n ( default " << GENERAL_DEFAULT_MIN_SAMPLES << " ) samples will be removed" << endl;
      cout << " -" << seed_s << " / --" << seed_l << setw( maxWidth - seed_l.size() ) 
	   << " " << "Seed (positive integer) for the Mersenne Twister random number generator" << endl;
      cout << endl;

      cout << "STOCHASTIC FOREST ARGUMENTS:" << endl;
      cout << " -" << nTrees_s << " / --" << nTrees_l << setw( maxWidth - nTrees_l.size() )
	   << " " << "Number of trees in the forest" << endl;
      cout << " -" << mTry_s << " / --" << mTry_l << setw( maxWidth - mTry_l.size() )
	   << " " << "Fraction of randomly drawn features per node split" << endl;
      cout << " -" << nMaxLeaves_s << " / --" << nMaxLeaves_l << setw( maxWidth - nMaxLeaves_l.size() )
	   << " " << "Maximum number of leaves per tree" << endl;
      cout << " -" << nodeSize_s << " / --" << nodeSize_l << setw( maxWidth - nodeSize_l.size() )
	   << " " << "Minimum number of train samples per node, affects tree depth" << endl;
      cout << " -" << shrinkage_s << " / --" << shrinkage_l << setw( maxWidth - shrinkage_l.size() )
	   << " " << "[GBT only] Shrinkage applied to evolving the residual" << endl;
      cout << endl;

      cout << "FILTER ARGUMENTS:" << endl;
      cout << " -" << isFilter_s << " / --" << isFilter_l << setw( maxWidth - isFilter_l.size() )
           << " " << "Set this flag if you want to perform feature selection (default OFF)" << endl;
      cout << " -" << nPerms_s << " / --" << nPerms_l << setw( maxWidth - nPerms_l.size() )
           << " " << "[Filter only] Number of permutations in statistical test (default " << ST_DEFAULT_N_PERMS << ")" << endl;
      cout << " -" << pValueThreshold_s << " / --" << pValueThreshold_l << setw( maxWidth - pValueThreshold_l.size() )
           << " " << "[Filter only] P-value threshold in statistical test (default " << ST_DEFAULT_P_VALUE_THRESHOLD << ")" << endl;
      cout << " -" << noSort_s << " / --" << noSort_l << setw( maxWidth - noSort_l.size() )
	   << " " << "[Filter only] Set this flag if you want the associations be listed unsorted (default OFF)" << endl; 
      cout << " -" << contrastOutput_s << " / --" << contrastOutput_l << setw( maxWidth - contrastOutput_l.size() )
	   << " " << "[Filter only] Output file for contrast value summary(ies) (one contrast value per permutation)" << endl;
      cout << endl;

      cout << "EXAMPLES:" << endl << endl;
      
      cout << "Performing feature selection using 'target' as the target variable:" << endl
	   << "bin/rf-ace --filter -I data.arff -i target -O associations.tsv" << endl << endl;

      cout << "Performing feature selection with 50 permutations and p-value threshold of 0.001:" << endl
	   << "bin/rf-ace --filter -I data.arff -i 5 -p 50 -t 0.001 -O associations.tsv" << endl << endl;

      cout << "Performing feature selection with regular Random Forest (p == 1) and listing associations unsorted:" << endl
	   << "bin/rf-ace --filter -I data.arff -i 5 -p 1 --noSort -O associations_unsorted.tsv" << endl << endl;
      
      cout << "Building Random Forest with 1000 trees and mTry of 10 and predicting with test data:" << endl
	   << "bin/rf-ace -I data.arff -i target -T testdata.arff -n 1000 -m 10 -O predictions.tsv" << endl << endl;

      cout << "Building Random Forest predictor and saving it to a file for later use:" << endl
	   << "bin/rf-ace -I data.arff -i target -O rf_predictor.sf" << endl << endl;

      cout << "Loading Random Forest predictor from file and predicting with test data:" << endl
	   << "bin/rf-ace -I predictor.sf -T testdata.arff -O predictions.tsv" << endl << endl;
      
    }
    
    void setRFDefaults() {

      nTrees                = RF_DEFAULT_N_TREES;
      mTry                  = RF_DEFAULT_M_TRY;
      nMaxLeaves            = RF_DEFAULT_N_MAX_LEAVES;
      nodeSize              = RF_DEFAULT_NODE_SIZE;
      shrinkage             = RF_DEFAULT_SHRINKAGE;

    }

    void setGBTDefaults() {

      nTrees                = GBT_DEFAULT_N_TREES;
      mTry                  = GBT_DEFAULT_M_TRY;
      nMaxLeaves            = GBT_DEFAULT_N_MAX_LEAVES;
      nodeSize              = GBT_DEFAULT_NODE_SIZE;
      shrinkage             = GBT_DEFAULT_SHRINKAGE;

    }
    
  };
  
}

#endif
