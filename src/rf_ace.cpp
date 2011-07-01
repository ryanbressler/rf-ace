#include <cstdlib>
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <stdio.h>

#include "getopt_pp.h"
#include "randomforest.hpp"
//#include "GBT.hpp"
#include "treedata.hpp"
#include "datadefs.hpp"

using namespace std;
using datadefs::num_t;
using namespace GetOpt;

enum MatrixFormat { FEATURE_ROWS, FEATURE_COLUMNS };

const int TARGETIDX_NOT_SET = -1;
const size_t DEFAULT_NTREES = 0;
const size_t DEFAULT_MTRY = 0;
const size_t DEFAULT_NODESIZE = 0;
const size_t DEFAULT_NPERMS = 50;
const num_t DEFAULT_PVALUETHRESHOLD = 0.10;

void printHeader()
{
  cout << endl;
  cout << " --------------------------------------------------------------- " << endl;
  cout << "| RF-ACE -- efficient feature selection with heterogeneous data |" << endl;
  cout << "|                                                               |" << endl;
  cout << "|  Version:      RF-ACE v0.4.0, July 1st, 2011                  |" << endl;
  cout << "|  Project page: http://code.google.com/p/rf-ace                |" << endl;
  cout << "|  Contact:      timo.p.erkkila@tut.fi                          |" << endl;
  cout << "|                                                               |" << endl;
  cout << "|              DEVELOPMENT VERSION, BUGS EXIST!                 |" << endl;
  cout << " --------------------------------------------------------------- " << endl;
}

void printHelp()
{
  cout << endl;
  cout << "REQUIRED ARGUMENTS:" << endl;
  cout << " -I / --input        input feature file (AFM or ARFF)" << endl;
  cout << " -i / --targetidx    target index, reference to order of appearance in the input file" << endl;
  cout << " -O / --output       output association file" << endl;
  cout << endl;
  cout << "OPTIONAL ARGUMENTS:" << endl;
  cout << " -n / --ntrees       number of trees per RF (default nsamples/nrealsamples)" << endl;
  cout << " -m / --mtry         number of randomly drawn features per node split (default sqrt(nfeatures))" << endl;
  cout << " -s / --nodesize     minimum number of train samples per node, affects tree depth (default max{5,nsamples/20})" << endl;
  cout << " -p / --nperms       number of Random Forests (default " << DEFAULT_NPERMS << ")" << endl;
  cout << " -t / --pthreshold   p-value threshold below which associations are listed (default " << DEFAULT_PVALUETHRESHOLD << ")" << endl;
  cout << endl;
}

void printHelpHint()
{
  cout << endl;
  cout << "To get started, type \"-h\" or \"--help\"" << endl;
}

int main(int argc, char* argv[])
{

  printHeader();

  if(argc == 1)
    {
      printHelp();
      return(EXIT_SUCCESS);
    }

  if(argc == 2)
    {
      string unaryHandle(argv[1]);
      if(unaryHandle == "-h" || unaryHandle == "--help")
	{
	  printHelp();
	  return(EXIT_SUCCESS);
	}
      else
	{
	  printHelpHint();
	  return(EXIT_FAILURE);
	}
    }

  string input = "";
  int targetIdx = TARGETIDX_NOT_SET;
  size_t nTrees = DEFAULT_NTREES;
  size_t mTry = DEFAULT_MTRY;
  size_t nodeSize = DEFAULT_NODESIZE;
  size_t nPerms = DEFAULT_NPERMS;
  num_t pValueThreshold = DEFAULT_PVALUETHRESHOLD;
  string output = "";

  GetOpt_pp ops(argc, argv);
        
  ops >> Option('I',"input",input);
  ops >> Option('i',"targetidx",targetIdx);
  ops >> Option('n',"ntrees",nTrees);
  ops >> Option('m',"mtry",mTry);
  ops >> Option('s',"nodesize",nodeSize);
  ops >> Option('p',"nperms",nPerms);
  ops >> Option('t',"pthreshold",pValueThreshold);
  ops >> Option('O',"output",output);

  if(input == "")
    {
      cerr << "Input file not specified" << endl;
      printHelpHint();
      return EXIT_FAILURE;
    }

  if(targetIdx == TARGETIDX_NOT_SET)
    {
      cerr << "targetidx (-i/--targetidx) must be set" << endl;
    }

  if(output == "")
    {
      cerr << "Output file not specified" << endl;
      printHelpHint();
      return EXIT_FAILURE;
    }

  //Read data into Treedata object
  cout << endl << "Reading file '" << input << "'" << endl;
  Treedata treedata(input);

  if(targetIdx < 0 || targetIdx >= static_cast<int>(treedata.nFeatures()))
    {
      cerr << "targetidx must point to a valid feature (0.." << treedata.nFeatures()-1 << ")" << endl;
      printHelpHint();
      return EXIT_FAILURE;
    }

  size_t nRealSamples = treedata.nRealSamples(targetIdx);

  if(nRealSamples == 0)
    {
      cout << "Target has no real samples, quitting..." << endl;
      return EXIT_SUCCESS;
    }

  num_t realFraction = 1.0*nRealSamples / treedata.nSamples();

  if(nTrees == DEFAULT_NTREES)
    {
      nTrees = static_cast<size_t>( 1.0*treedata.nSamples() / realFraction );
    }
  
  if(mTry == DEFAULT_MTRY)
    {
      mTry = static_cast<size_t>( floor( sqrt( 1.0*treedata.nFeatures() ) ) );   
    }

  if(nodeSize == DEFAULT_NODESIZE)
    {
      nodeSize = 5;
      size_t altNodeSize = static_cast<size_t>( ceil( 1.0*nRealSamples/20 ) );
      if(altNodeSize > nodeSize)
	{
	  nodeSize = altNodeSize;
	}
    }
  
  cout << endl;
  cout << "RF-ACE parameter configuration:" << endl;
  cout << "  --input      = " << input << endl;
  cout << "  --nsamples   = " << nRealSamples << " / " << treedata.nSamples() << " (" << 100 * ( 1 - realFraction )<< "% missing)" << endl;
  cout << "  --nfeatures  = " << treedata.nFeatures() << endl;
  cout << "  --targetidx  = " << targetIdx << ", header '" << treedata.getFeatureName(targetIdx) << "'" << endl;
  cout << "  --ntrees     = " << nTrees << endl;
  cout << "  --mtry       = " << mTry << endl;
  cout << "  --nodesize   = " << nodeSize << endl;
  cout << "  --nperms     = " << nPerms << endl;
  cout << "  --pthresold  = " << pValueThreshold << endl;
  //cout << "  --alpha      = " << alpha << endl;
  cout << "  --output     = " << output << endl << endl;

  assert(treedata.nFeatures() >= mTry);
  assert(treedata.nSamples() > 2*nodeSize);

  Randomforest RF(&treedata,targetIdx,nTrees,mTry,nodeSize);
  //RF.setTarget(targetIdx);
  //treedata.print();
  //RF.blacklist_and_kill(0.8,blistheaders,blistcorrelations);

  //size_t nperms = 9;
  //num_t alpha = 0.50;
  vector<num_t> pValues(treedata.nFeatures());
  vector<num_t> importanceValues(treedata.nFeatures());
  //vector<num_t> ivalues(treedata.nfeatures());
  
  cout << "Growing " << nPerms << " Random Forests (RFs), please wait..." << endl;
  RF.learn(nPerms,pValues,importanceValues);
  //cout << "Time elapsed: " << float(clock() - time_start)/CLOCKS_PER_SEC << " seconds" << endl;
  
  vector<size_t> refIcs(treedata.nFeatures());
  //vector<string> fnames = treedata.featureheaders();
  datadefs::sortDataAndMakeRef(pValues,refIcs);
  datadefs::sortFromRef<num_t>(importanceValues,refIcs);
  //targetIdx = refIcs[targetIdx];
  
  string targetName = treedata.getFeatureName(targetIdx);

  if(pValues[0] <= pValueThreshold)
    {
      //cout << "Writing associations to file '" << output << "'" << endl;
      //ofstream os(output.c_str());
      FILE* po;
      po = fopen(output.c_str(),"w");
      
      for(size_t featureIdx = 0; featureIdx < treedata.nFeatures(); ++featureIdx)
	{

	  if(pValues[featureIdx] > pValueThreshold)
	    {
	      break;
	    }

	  if(refIcs[featureIdx] == static_cast<size_t>(targetIdx))
	    {
	      //cout << refIcs[featureIdx] << " == " << targetIdx << " (" << targetHeader << ")" << endl;
	      continue;
	    }


	  fprintf(po,"%s\t%s\t%9.8f\t%9.8f\t%9.8f\n",targetName.c_str(),treedata.getFeatureName(refIcs[featureIdx]).c_str(),pValues[featureIdx],importanceValues[featureIdx],treedata.pearsonCorrelation(targetIdx,refIcs[featureIdx]));

	  //os << target_str << "\t" << treedata.get_featureheader(ref_ics[i]) << "\t" 
	  //   << pvalues[i] << "\t" << ivalues[i] << "\t" << treedata.corr(targetidx,ref_ics[i]) << endl;
	}
      //os.close();
      fclose(po);
      cout << endl << "Association file created. Format:" << endl;
      cout << "TARGET   PREDICTOR   P-VALUE   IMPORTANCE   CORRELATION" << endl << endl;
      cout << "Done." << endl;
    }
  else
    {
      cout << endl << "No significant associations found, quitting..." << endl;
    }

  
 
  return(EXIT_SUCCESS);
}
