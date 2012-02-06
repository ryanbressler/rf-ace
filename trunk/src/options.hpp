#include <cstdlib>
#include <string>
#include <sstream>
#include <iostream>
#include "argparse.hpp"
#include "datadefs.hpp"

using namespace std;
using datadefs::num_t;

namespace options {

  void printHeader(ostream& output) {

    output << endl;
    output << " ------------------------------------------------------- " << endl;
    output << "|  RF-ACE version:  0.9.9, February 2nd, 2012           |" << endl;
    output << "|    Project page:  http://code.google.com/p/rf-ace     |" << endl;
    output << "|     Report bugs:  timo.p.erkkila@tut.fi               |" << endl;
    output << " ------------------------------------------------------- " << endl;
    output << endl;
  }

  void printHelpHint() {
    cout << endl;
    cout << "To get started, type \"-h\" or \"--help\"" << endl;
  }

  const bool   GENERAL_DEFAULT_PRINT_HELP = false;
  const bool   GENERAL_DEFAULT_NO_PREDICTION = false; // TEMPORARY VARIABLE
  const num_t  GENERAL_DEFAULT_P_VALUE_THRESHOLD = 0.05;
  const bool   GENERAL_DEFAULT_NO_FILTER = false;
  const char   GENERAL_DEFAULT_DATA_DELIMITER = '\t';
  const char   GENERAL_DEFAULT_HEADER_DELIMITER = ':';
  const size_t GENERAL_DEFAULT_MIN_SAMPLES = 5;
  
  const size_t RF_DEFAULT_N_TREES = 1000; // zero means it will be estimated from the data by default
  const size_t RF_DEFAULT_M_TRY = 0; // same here ...
  const size_t RF_DEFAULT_N_MAX_LEAVES = 100;
  const size_t RF_DEFAULT_NODE_SIZE = 3; // ... and here
  const size_t RF_DEFAULT_N_PERMS = 20;
  const num_t  RF_DEFAULT_P_VALUE_THRESHOLD = 0.05;
  
  const size_t GBT_DEFAULT_N_TREES = 100;
  const size_t GBT_DEFAULT_N_MAX_LEAVES = 6;
  const num_t  GBT_DEFAULT_SHRINKAGE = 0.1;
  const num_t  GBT_DEFAULT_SUB_SAMPLE_SIZE = 0.5;

  struct Filte_options {

          size_t nTrees;
    const string nTrees_s;
    const string nTrees_l;

          size_t mTry;
    const string mTry_s;
    const string mTry_l;

          size_t nMaxLeaves;
    const string nMaxLeaves_s;
    const string nMaxLeaves_l;

          size_t nodeSize;
    const string nodeSize_s;
    const string nodeSize_l;

          size_t nPerms;
    const string nPerms_s;
    const string nPerms_l;

    num_t pValueThreshold;
    const string pValueThreshold_s;
    const string pValueThreshold_l;

    Filte_options(const int argc, char* const argv[]):
    
      nTrees(RF_DEFAULT_N_TREES),
      nTrees_s("n"),
      nTrees_l("RF_ntrees"),

      mTry(RF_DEFAULT_M_TRY),
      mTry_s("m"),
      mTry_l("RF_mtry"),

      nMaxLeaves(RF_DEFAULT_N_MAX_LEAVES),
      nMaxLeaves_s("a"),
      nMaxLeaves_l("RF_maxleaves"),

      nodeSize(RF_DEFAULT_NODE_SIZE),
      nodeSize_s("s"),
      nodeSize_l("RF_nodesize"),

      nPerms(RF_DEFAULT_N_PERMS),
      nPerms_s("p"),
      nPerms_l("RF_nperms"),

      pValueThreshold(RF_DEFAULT_P_VALUE_THRESHOLD),
      pValueThreshold_s("t"),
      pValueThreshold_l("pthreshold") {

      // Read the user parameters ...
      ArgParse parser(argc,argv);

      parser.getArgument<size_t>(nTrees_s,nTrees_l,nTrees);
      parser.getArgument<size_t>(mTry_s, mTry_l, mTry);
      parser.getArgument<size_t>(nMaxLeaves_s, nMaxLeaves_l, nMaxLeaves);
      parser.getArgument<size_t>(nodeSize_s, nodeSize_l, nodeSize);
      parser.getArgument<size_t>(nPerms_s, nPerms_l, nPerms);
      parser.getArgument(pValueThreshold_s, pValueThreshold_l, pValueThreshold);

    }

  };

  struct General_options {
            
            bool printHelp;
    const string printHelp_s;
    const string printHelp_l;
    
          string trainInput;
    const string trainInput_s;
    const string trainInput_l;
    
          string associationOutput;
    const string associationOutput_s;
    const string associationOutput_l;
    
          string targetStr;
    const string targetStr_s;
    const string targetStr_l;
    
          string whiteListInput;
    const string whiteListInput_s;
    const string whiteListInput_l;

          string blackListInput;
    const string blackListInput_s;
    const string blackListInput_l;
    
          string testInput;
    const string testInput_s;
    const string testInput_l;
    
          string predictionOutput;
    const string predictionOutput_s;
    const string predictionOutput_l;
    
          string logOutput;
    const string logOutput_s;
    const string logOutput_l;
    
          string forestOutput;
    const string forestOutput_s;
    const string forestOutput_l;
    
           num_t pValueThreshold;
    const string pValueThreshold_s;
    const string pValueThreshold_l;
    
            char dataDelimiter;
    const string dataDelimiter_s;
    const string dataDelimiter_l;
    
            char headerDelimiter;
    const string headerDelimiter_s;
    const string headerDelimiter_l;

          size_t pruneFeatures;
    const string pruneFeatures_s;
    const string pruneFeatures_l;
    
            bool noFilter;
    const string noFilter_s;
    const string noFilter_l;
    
    General_options(const int argc, char* const argv[]):
      printHelp(GENERAL_DEFAULT_PRINT_HELP),
      printHelp_s("h"),
      printHelp_l("help"),    
      
      trainInput(""),
      trainInput_s("I"),
      trainInput_l("traindata"),
      
      associationOutput(""),
      associationOutput_s("O"),
      associationOutput_l("associations"),
      
      targetStr(""),
      targetStr_s("i"),
      targetStr_l("target"),
      
      whiteListInput(""),
      whiteListInput_s("W"),
      whiteListInput_l("whitelist"),
      
      blackListInput(""),
      blackListInput_s("B"),
      blackListInput_l("blacklist"),

      testInput(""),
      testInput_s("T"),
      testInput_l("testdata"),
      
      predictionOutput(""),
      predictionOutput_s("P"),
      predictionOutput_l("predictions"),
      
      logOutput(""),
      logOutput_s("L"),
      logOutput_l("log"),
      
      forestOutput(""),
      forestOutput_s("F"),
      forestOutput_l("forest"),
      
      pValueThreshold(GENERAL_DEFAULT_P_VALUE_THRESHOLD),
      pValueThreshold_s("t"),
      pValueThreshold_l("pthreshold"),
      
      dataDelimiter(GENERAL_DEFAULT_DATA_DELIMITER),
      dataDelimiter_s("D"),
      dataDelimiter_l("data_delim"),
      
      headerDelimiter(GENERAL_DEFAULT_HEADER_DELIMITER),
      headerDelimiter_s("H"),
      headerDelimiter_l("head_delim"),
      
      pruneFeatures(GENERAL_DEFAULT_MIN_SAMPLES),
      pruneFeatures_s("X"),
      pruneFeatures_l("prune_features"),

      noFilter(GENERAL_DEFAULT_NO_FILTER),
      noFilter_s("N"),
      noFilter_l("no_filter") {

      // Read the user parameters ...
      ArgParse parser(argc,argv);

      // First read general options
      parser.getFlag(printHelp_s, printHelp_l, printHelp);
      parser.getArgument<string>(trainInput_s, trainInput_l, trainInput);
      parser.getArgument<string>(targetStr_s, targetStr_l, targetStr);
      parser.getArgument<string>(associationOutput_s, associationOutput_l, associationOutput);
      parser.getArgument<string>(whiteListInput_s, whiteListInput_l, whiteListInput);
      parser.getArgument<string>(blackListInput_s, blackListInput_l, blackListInput);
      parser.getArgument<string>(testInput_s, testInput_l, testInput);
      parser.getArgument<string>(predictionOutput_s, predictionOutput_l, predictionOutput);
      parser.getArgument<string>(logOutput_s,logOutput_l,logOutput);
      parser.getArgument<string>(forestOutput_s,forestOutput_l,forestOutput);
      parser.getArgument<num_t>(pValueThreshold_s, pValueThreshold_l, pValueThreshold);
      string dataDelimiter,headerDelimiter;
      parser.getArgument<string>(dataDelimiter_s, dataDelimiter_l, dataDelimiter);
      parser.getArgument<string>(headerDelimiter_s, headerDelimiter_l, headerDelimiter);
      parser.getArgument<size_t>(pruneFeatures_s, pruneFeatures_l, pruneFeatures);
      parser.getFlag(noFilter_s, noFilter_l, noFilter);
      stringstream ss(dataDelimiter);
      ss >> dataDelimiter;
      //cout << dataDelimiter;

      ss.clear();
      ss.str("");
      ss << headerDelimiter;
      ss >> headerDelimiter;
 
    }
  };
  
  struct RF_options {
    
          size_t nTrees;
    const string nTrees_s;
    const string nTrees_l;
    
          size_t mTry;
    const string mTry_s;
    const string mTry_l;
    
          size_t nMaxLeaves;
    const string nMaxLeaves_s;
    const string nMaxLeaves_l;
    
          size_t nodeSize;
    const string nodeSize_s;
    const string nodeSize_l;
    
          size_t nPerms;
    const string nPerms_s;
    const string nPerms_l;

    num_t pValueThreshold;
    const string pValueThreshold_s;
    const string pValueThreshold_l;
    
    RF_options(const int argc, char* const argv[]):
      nTrees(RF_DEFAULT_N_TREES),
      nTrees_s("n"),
      nTrees_l("RF_ntrees"),
      
      mTry(RF_DEFAULT_M_TRY),
      mTry_s("m"),
      mTry_l("RF_mtry"),
      
      nMaxLeaves(RF_DEFAULT_N_MAX_LEAVES),
      nMaxLeaves_s("a"),
      nMaxLeaves_l("RF_maxleaves"),
      
      nodeSize(RF_DEFAULT_NODE_SIZE),
      nodeSize_s("s"),
      nodeSize_l("RF_nodesize"),
      
      nPerms(RF_DEFAULT_N_PERMS),
      nPerms_s("p"),
      nPerms_l("RF_nperms"),

      pValueThreshold(RF_DEFAULT_P_VALUE_THRESHOLD),
      pValueThreshold_s("t"),
      pValueThreshold_l("pthreshold") {
      
      // Read the user parameters ...
      ArgParse parser(argc,argv);

      parser.getArgument<size_t>(nTrees_s,nTrees_l,nTrees);
      parser.getArgument<size_t>(mTry_s, mTry_l, mTry);
      parser.getArgument<size_t>(nMaxLeaves_s, nMaxLeaves_l, nMaxLeaves);
      parser.getArgument<size_t>(nodeSize_s, nodeSize_l, nodeSize);
      parser.getArgument<size_t>(nPerms_s, nPerms_l, nPerms);
      parser.getArgument(pValueThreshold_s, pValueThreshold_l, pValueThreshold);

    }
  };
  
  struct GBT_options {
    
          size_t nTrees;
    const string nTrees_s;
    const string nTrees_l;
    
          size_t nMaxLeaves;
    const string nMaxLeaves_s;
    const string nMaxLeaves_l;
    
           num_t shrinkage;
    const string shrinkage_s;
    const string shrinkage_l;
    
           num_t subSampleSize;
    const string subSampleSize_s;
    const string subSampleSize_l;
    
    GBT_options(const int argc, char* const argv[]):    
      nTrees(GBT_DEFAULT_N_TREES),
      nTrees_s("r"),
      nTrees_l("GBT_ntrees"),
      
      nMaxLeaves(GBT_DEFAULT_N_MAX_LEAVES),
      nMaxLeaves_s("l"),
      nMaxLeaves_l("GBT_maxleaves"),
      
      shrinkage(GBT_DEFAULT_SHRINKAGE),
      shrinkage_s("z"),
      shrinkage_l("GBT_shrinkage"),
      
      subSampleSize(GBT_DEFAULT_SUB_SAMPLE_SIZE),
      subSampleSize_s("u"),
      subSampleSize_l("GBT_samplesize") {

      // Read the user parameters ...
      ArgParse parser(argc,argv);      
      
      parser.getArgument<size_t>(nTrees_s, nTrees_l, nTrees);
      parser.getArgument<size_t>(nMaxLeaves_s, nMaxLeaves_l, nMaxLeaves);
      parser.getArgument<num_t>(shrinkage_s, shrinkage_l, shrinkage);
      parser.getArgument<num_t>(subSampleSize_s, subSampleSize_l, subSampleSize);
      
    }
  };

  void printHelp(const General_options& geno, const RF_options& rfo, const GBT_options& gbto) {
    
    size_t maxwidth = 20;
    cout << endl;
    
    cout << "REQUIRED ARGUMENTS:" << endl;
    cout << " -" << geno.trainInput_s << " / --" << geno.trainInput_l << setw( maxwidth - geno.trainInput_l.size() )
	 << " " << "Train data input file, associations will be sought from this data. Supported formats: AFM and ARFF" << endl;
    cout << " -" << geno.targetStr_s << " / --" << geno.targetStr_l << setw( maxwidth - geno.targetStr_l.size() )
	 << " " << "Target, specified as integer or string that is to be matched with the content of input" << endl;
    cout << endl;

    cout << "OPTIONAL ARGUMENTS:" << endl;
    cout << " -" << geno.associationOutput_s << " / --" << geno.associationOutput_l << setw( maxwidth - geno.associationOutput_l.size() )
	 << " " << "Association output file" << endl;
    cout << " -" << geno.testInput_s << " / --" << geno.testInput_l << setw( maxwidth - geno.testInput_l.size() )
	 << " " << "Test data input file, predictions will be made from this data" << endl;
    cout << " -" << geno.predictionOutput_s << " / --" << geno.predictionOutput_l << setw( maxwidth - geno.predictionOutput_l.size() )
	 << " " << "Prediction output file" << endl;
    cout << " -" << geno.logOutput_s << " / --" << geno.logOutput_l << setw( maxwidth - geno.logOutput_l.size() )
	 << " " << "Log output file" << endl;
    cout << " -" << geno.forestOutput_s << " / --" << geno.forestOutput_l << setw( maxwidth - geno.forestOutput_l.size() )
	 << " " << "Forest output file" << endl;
    cout << " -" << geno.whiteListInput_s << " / --" << geno.whiteListInput_l << setw( maxwidth - geno.whiteListInput_l.size() )
	 << " " << "White list of features to be included in the input file(s)." << endl;
    cout << " -" << geno.blackListInput_s << " / --" << geno.blackListInput_l << setw( maxwidth - geno.blackListInput_l.size() )
         << " " << "Black list of features to be excluded from the input file(s)." << endl;
    cout << " -" << geno.dataDelimiter_s << " / --" << geno.dataDelimiter_l << setw( maxwidth - geno.dataDelimiter_l.size() )
	 << " " << "Data delimiter (default \\t)" << endl;
    cout << " -" << geno.headerDelimiter_s << " / --" << geno.headerDelimiter_l << setw( maxwidth - geno.headerDelimiter_l.size() )
	 << " " << "Header delimiter that separates the N and C symbols in feature headers from the rest (default " << GENERAL_DEFAULT_HEADER_DELIMITER << ")" << endl;
    cout << " -" << geno.pruneFeatures_s << " / --" << geno.pruneFeatures_l << setw( maxwidth - geno.pruneFeatures_l.size() )
	 << " " << "Features with less than n ( default " << GENERAL_DEFAULT_MIN_SAMPLES << " ) samples will be removed" << endl;
    cout << " -" << geno.noFilter_s << " / --" << geno.noFilter_l << setw( maxwidth - geno.noFilter_l.size() )
	 << " " << "Set this flag if you want to omit applying feature filtering with RFs" << endl;
    cout << endl;

    cout << "OPTIONAL ARGUMENTS -- RANDOM FOREST:" << endl;
    cout << " -" << rfo.nTrees_s << " / --" << rfo.nTrees_l << setw( maxwidth - rfo.nTrees_l.size() )
	 << " " << "Number of trees per RF (default " << RF_DEFAULT_N_TREES << ")" << endl;
    cout << " -" << rfo.mTry_s << " / --" << rfo.mTry_l << setw( maxwidth - rfo.mTry_l.size() )
	 << " " << "Number of randomly drawn features per node split (default floor(0.1*nFeatures))" << endl;
    cout << " -" << rfo.nMaxLeaves_s << " / --" << rfo.nMaxLeaves_l << setw( maxwidth - rfo.nMaxLeaves_l.size() )
	 << " " << "Maximum number of leaves per tree (default " << RF_DEFAULT_N_MAX_LEAVES << ")" << endl;
    cout << " -" << rfo.nodeSize_s << " / --" << rfo.nodeSize_l << setw( maxwidth - rfo.nodeSize_l.size() )
	 << " " << "Minimum number of train samples per node, affects tree depth (default " << RF_DEFAULT_NODE_SIZE << ")" << endl;
    cout << " -" << rfo.nPerms_s << " / --" << rfo.nPerms_l << setw( maxwidth - rfo.nPerms_l.size() )
	 << " " << "Number of Random Forests (default " << RF_DEFAULT_N_PERMS << ")" << endl;
    cout << endl;

    cout << "OPTIONAL ARGUMENTS -- FEATURE FILTER:" << endl;
    cout << " -" << geno.pValueThreshold_s << " / --" << geno.pValueThreshold_l << setw( maxwidth - geno.pValueThreshold_l.size() )
	 << " " << "p-value threshold below which associations are listed (default "
	 << GENERAL_DEFAULT_P_VALUE_THRESHOLD << ")" << endl;
    cout << endl;

    cout << "OPTIONAL ARGUMENTS -- GRADIENT BOOSTING TREES:" << endl;
    cout << " -" << gbto.nTrees_s << " / --" << gbto.nTrees_l << setw( maxwidth - gbto.nTrees_l.size() )
	 << " " << "Number of trees in the GBT (default " << GBT_DEFAULT_N_TREES << ")" << endl;
    cout << " -" << gbto.nMaxLeaves_s << " / --" << gbto.nMaxLeaves_l << setw( maxwidth - gbto.nMaxLeaves_l.size() )
	 << " " << "Maximum number of leaves per tree (default " << GBT_DEFAULT_N_MAX_LEAVES << ")" << endl;
    cout << " -" << gbto.shrinkage_s << " / --" << gbto.shrinkage_l << setw( maxwidth - gbto.shrinkage_l.size() )
	 << " " << "Shrinkage applied to evolving the residual (default " << GBT_DEFAULT_SHRINKAGE << ")" << endl;
    cout << " -" << gbto.subSampleSize_s << " / --" << gbto.subSampleSize_l << setw( maxwidth - gbto.subSampleSize_l.size() )
	 << " " << "Sample size fraction for training the trees (default " << GBT_DEFAULT_SUB_SAMPLE_SIZE << ")" << endl;
    cout << endl;

    cout << "EXAMPLES:" << endl;
    cout << endl;

    cout << "List features associated with feature called \"C:AFFECTION\":" << endl;
    cout << "bin/rf_ace --traindata featurematrix.afm --target C:AFFECTION --associations associations.tsv" << endl << endl;

    cout << "Train the model for \"C:AFFECTION\" with \"original.afm\", and predict new data for \"C:AFFECTION\" with \"newdata.afm\":" << endl;
    cout << "bin/rf_ace --traindata original.afm --target C:AFFECTION --testdata newdata.afm --predictions.tsv" << endl << endl;
  }
 
}
