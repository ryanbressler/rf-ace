#include <cstdlib>
#include <string>
#include "datadefs.hpp"

using namespace std;
using datadefs::num_t;

namespace options {
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
  
  const size_t GBT_DEFAULT_N_TREES = 100;
  const size_t GBT_DEFAULT_N_MAX_LEAVES = 6;
  const num_t  GBT_DEFAULT_SHRINKAGE = 0.1;
  const num_t  GBT_DEFAULT_SUB_SAMPLE_SIZE = 0.5;
  
  struct General_options {
    
    bool   printHelp;
    string printHelp_s;
    string printHelp_l;
    
    string trainInput;
    string trainInput_s;
    string trainInput_l;
    
    string associationOutput;
    string associationOutput_s;
    string associationOutput_l;
    
    string targetStr;
    string targetStr_s;
    string targetStr_l;
    
    string whiteListInput;
    string whiteListInput_s;
    string whiteListInput_l;

    string blackListInput;
    string blackListInput_s;
    string blackListInput_l;
    
    string testInput;
    string testInput_s;
    string testInput_l;
    
    string predictionOutput;
    string predictionOutput_s;
    string predictionOutput_l;
    
    string logOutput;
    string logOutput_s;
    string logOutput_l;
    
    string forestOutput;
    string forestOutput_s;
    string forestOutput_l;
    
    num_t  pValueThreshold;
    string pValueThreshold_s;
    string pValueThreshold_l;
    
    char   dataDelimiter;
    string dataDelimiter_s;
    string dataDelimiter_l;
    
    char   headerDelimiter;
    string headerDelimiter_s;
    string headerDelimiter_l;

    size_t pruneFeatures;
    string pruneFeatures_s;
    string pruneFeatures_l;
    
    bool   noFilter;
    string noFilter_s;
    string noFilter_l;
    
    General_options():
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
      noFilter_l("no_filter") {}
  };
  
  struct RF_options {
    
    size_t nTrees;
    string nTrees_s;
    string nTrees_l;
    
    size_t mTry;
    string mTry_s;
    string mTry_l;
    
    size_t nMaxLeaves;
    string nMaxLeaves_s;
    string nMaxLeaves_l;
    
    size_t nodeSize;
    string nodeSize_s;
    string nodeSize_l;
    
    size_t nPerms;
    string nPerms_s;
    string nPerms_l;
    
    RF_options():
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
      nPerms_l("RF_nperms") {}
  };
  
  struct GBT_options {
    
    size_t nTrees;
    string nTrees_s;
    string nTrees_l;
    
    size_t nMaxLeaves;
    string nMaxLeaves_s;
    string nMaxLeaves_l;
    
    num_t shrinkage;
    string shrinkage_s;
    string shrinkage_l;
    
    num_t subSampleSize;
    string subSampleSize_s;
    string subSampleSize_l;
    
    GBT_options():    
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
      subSampleSize_l("GBT_samplesize") {}
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
