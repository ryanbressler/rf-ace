#include <cstdlib>
#include <iostream>
#include "argparse.hpp"
#include "stochasticforest.hpp"
#include "treedata.hpp"

using namespace std;


const size_t DEFAULT_TARGETIDX = 0;
const size_t DEFAULT_NTREES = 500;
const size_t DEFAULT_NODESIZE = 5;
const num_t DEFAULT_SHRINKAGE = 0.2;
const num_t DEFAULT_SUBSAMPLE = 0.5;

int main(const int argc, char* const argv[]) {

  //------------------------------------------------------------------------
  // 0: parameters
  if(argc == 1 || argc == 2) {
    if(argc == 2) {
      string helphandle(argv[1]);
      if (helphandle != "-h" && helphandle != "--help") {
        cerr << "use -h or --help to get started" << endl;
        return EXIT_FAILURE;
      }
    }

    cout << endl;
    cout << "REQUIRED ARGUMENTS:" << endl;
    cout << "-I / --input        input feature matrix" << endl;
    cout << "-O / --output       output association file" << endl;
    cout << endl;
    cout << "OPTIONAL ARGUMENTS:" << endl;
    cout << "-i / --targetidx    target index, ref. to feature matrix (default " << DEFAULT_TARGETIDX << ")" << endl;
    cout << "-n / --ntrees       number of trees per GBT forest (default " << DEFAULT_NTREES << ")" << endl;
    cout << "-s / --nodesize     minimum number of train samples per node, affects tree depth (default " << DEFAULT_NODESIZE << ")" << endl;
    cout << "-z / --shrinkage    shrinkage (default " << DEFAULT_SHRINKAGE << ")" << endl;
    cout << "-u / --subsample    subsample size (default " << DEFAULT_SUBSAMPLE << ")" << endl;
    cout << endl;
    return EXIT_SUCCESS;
  }

  cout << endl;
  cout << "  ----------------------------------" << endl;
  cout << "  ---  GBT_benchmark version 0.0.2    ---" << endl;
  cout << "  ----------------------------------" << endl;

  //using namespace GetOpt;
  string input = "";
  size_t targetIdx = DEFAULT_TARGETIDX;
  size_t ntrees = DEFAULT_NTREES;
  size_t nodesize = DEFAULT_NODESIZE;
  num_t shrinkage = DEFAULT_SHRINKAGE;
  num_t subSampleSize = DEFAULT_SUBSAMPLE;
  string output = "";

  ArgParse parser(argc,argv);
  parser.getArgument<string>("I","input",input);
  parser.getArgument<size_t>("i","target",targetIdx);
  parser.getArgument<size_t>("n","ntrees",ntrees);
  parser.getArgument<string>("O","output",output);
  parser.getArgument<num_t>("z","shrinkage",shrinkage);
  parser.getArgument<num_t>("u","subsample",subSampleSize);


  if(input == "") {
    cerr << "input file not specified" << endl;
    return EXIT_FAILURE;
  }

  if(output == "") {
    cerr << "output file not specified" << endl;
    return EXIT_FAILURE;
  }


  //------------------------------------------------------------------------
  // 1: read data into Treedata class (features are rows)
  cout <<endl<< "READ:" << endl;
  Treedata treeData(input);
  size_t nSamples = treeData.nSamples();


  //------------------------------------------------------------------------
  // 2: construct a GBT object
  cout <<endl<< "CONSTRUCT:" << endl;
  size_t nTrees(ntrees);
  size_t nMaxLeaves(nodesize);
  StochasticForest myGbtForest( &treeData, targetIdx, nTrees);


  //------------------------------------------------------------------------
  // 3: grow the forest
  cout <<endl<< "GROWING:" << endl;
  myGbtForest.learnGBT(nMaxLeaves, shrinkage, subSampleSize);

  //------------------------------------------------------------------------
  // 4: predict using the forest
  cout << "PREDICTION: "<<nSamples<<" samples. Target="<<targetIdx<<endl;
  vector<num_t>  confidence(nSamples);
  vector<string> prediction(nSamples);
  vector<num_t>  numPrediction(nSamples);
  if ( treeData.isFeatureNumerical(targetIdx) )
    myGbtForest.predict(numPrediction, confidence);
  else
    myGbtForest.predict(prediction, confidence);

  vector<num_t> target(nSamples);
  treeData.getFeatureData(targetIdx, target);

  // diagnostic print out the true and the prediction (on train data)
  if ( treeData.isFeatureNumerical(targetIdx) ) {
    for (size_t i=0; i<nSamples; i++) {
      cout << i << "\t" << target[i] << "\t" << numPrediction[i] <<"\t" << confidence[i]<<endl;
    }

  } else {

    size_t errors = 0;
    for (size_t i=0; i<nSamples; i++) {
      string trueStr = treeData.getRawFeatureData(targetIdx,i);
      if (trueStr != prediction[i]) ++errors;
      cout << i << "\t" << trueStr << "\t" << prediction[i] <<"\t" << confidence[i]<<endl;
    }
    cout << "Errors " <<errors<<"/"<<nSamples<<" = " << 100.0*errors/nSamples<<"%."<<endl;
  }

  return(EXIT_SUCCESS);
}
