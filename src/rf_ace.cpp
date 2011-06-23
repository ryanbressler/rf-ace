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

const size_t DEFAULT_TARGETIDX = 0;
const size_t DEFAULT_NTREES = 0;
const size_t DEFAULT_MTRY = 0;
const size_t DEFAULT_NODESIZE = 5;
const size_t DEFAULT_NPERMS = 20;
const num_t DEFAULT_PVALUETHRESHOLD = 0.10;
//const num_t DEFAULT_ALPHA = 0.95;

int main(int argc, char* argv[])
{

  if(argc == 1 || argc == 2)
    {

      if(argc == 2)
	{
	  string helpHandle(argv[1]);
	  if(helpHandle != "-h" && helpHandle != "--help")
	    {
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
      cout << "-n / --ntrees       number of trees per RF (default 2*nsamples/nrealsamples)" << endl;
      cout << "-m / --mtry         number of randomly drawn features per node split (default sqrt(nfeatures))" << endl;
      cout << "-s / --nodesize     minimum number of train samples per node, affects tree depth (default " << DEFAULT_NODESIZE << ")" << endl;
      cout << "-p / --nperms       number of Random Forests (default " << DEFAULT_NPERMS << ")" << endl;
      cout << "-t / --pthreshold   p-value threshold below which associations are listed (default " << DEFAULT_PVALUETHRESHOLD << ")" << endl;
      //cout << "-a / --alpha        percentile of contrast importances, defines stringency of the t-test (default " << DEFAULT_ALPHA << ")" << endl;
      cout << endl;
      return EXIT_SUCCESS;
    }

  cout << endl;
  cout << "  ----------------------------------" << endl;
  cout << "  ---    RF-ACE version 0.3.0    ---" << endl;
  cout << "  ----------------------------------" << endl;

  //using namespace GetOpt;
  string input = "";
  size_t targetIdx = DEFAULT_TARGETIDX;
  size_t nTrees = DEFAULT_NTREES;
  size_t mTry = DEFAULT_MTRY;
  size_t nodeSize = DEFAULT_NODESIZE;
  size_t nPerms = DEFAULT_NPERMS;
  num_t pValueThreshold = DEFAULT_PVALUETHRESHOLD;
  //num_t alpha = DEFAULT_ALPHA;
  //bool is_featurerows = true;
  string output = "";

  GetOpt_pp ops(argc, argv);
        
  ops >> Option('I',"input",input);
  ops >> Option('i',"targetidx",targetIdx);
  ops >> Option('n',"ntrees",nTrees);
  ops >> Option('m',"mtry",mTry);
  ops >> Option('p',"nperms",nPerms);
  ops >> Option('t',"pthreshold",pValueThreshold);
  ops >> Option('O',"output",output);

  if(input == "")
    {
      cerr << "input file not specified" << endl;
      return EXIT_FAILURE;
    }

  if(output == "")
    {
      cerr << "output file not specified" << endl;
      return EXIT_FAILURE;
    }

  //Read data into Treedata object
  cout << endl << "Reading file '" << input << "'" << endl;
  Treedata treedata(input);

 
  if(targetIdx >= treedata.nfeatures())
    {
      cerr << "targetidx must point to a valid feature (0.." << treedata.nfeatures()-1 << ")" << endl;
      return EXIT_FAILURE;
    }

  size_t nRealSamples = treedata.nrealsamples(targetIdx);

  if(nRealSamples == 0)
    {
      cout << "Target has no real samples, quitting..." << endl;
      return EXIT_SUCCESS;
    }

  num_t realFraction = static_cast<num_t>(nRealSamples) / static_cast<num_t>(treedata.nsamples());

  if(nTrees == DEFAULT_NTREES)
    {
      nTrees = static_cast<size_t>(2.0 * static_cast<num_t>(treedata.nsamples()) / realFraction );
    }
  
  if(mTry == DEFAULT_MTRY)
    {
      mTry = static_cast<size_t>(floor(sqrt(static_cast<num_t>(treedata.nfeatures()))));   
    }

  //treedata.print(targetidx);
  
  cout << endl;
  cout << "RF-ACE parameter configuration:" << endl;
  cout << "  --input      = " << input << endl;
  cout << "  --nsamples   = " << nRealSamples << " / " << treedata.nsamples() << " (" << 100 * ( 1 - realFraction )<< "% missing)" << endl;
  cout << "  --nfeatures  = " << treedata.nfeatures() << endl;
  cout << "  --targetidx  = " << targetIdx << ", header '" << treedata.get_featureheader(targetIdx) << "'" << endl;
  cout << "  --ntrees     = " << nTrees << endl;
  cout << "  --mtry       = " << mTry << endl;
  cout << "  --nodesize   = " << nodeSize << endl;
  cout << "  --nperms     = " << nPerms << endl;
  cout << "  --pthresold  = " << pValueThreshold << endl;
  //cout << "  --alpha      = " << alpha << endl;
  cout << "  --output     = " << output << endl << endl;

  assert(treedata.nfeatures() >= mTry);
  assert(treedata.nsamples() > 2*nodeSize);

  Randomforest RF(&treedata,nTrees,mTry,nodeSize);
  RF.select_target(targetIdx);
  //RF.blacklist_and_kill(0.8,blistheaders,blistcorrelations);

  //size_t nperms = 9;
  //num_t alpha = 0.50;
  vector<num_t> pValues(treedata.nfeatures());
  vector<num_t> importanceValues(treedata.nfeatures());
  //vector<num_t> ivalues(treedata.nfeatures());
  
  cout << "Growing " << nPerms << " Random Forests (RFs), please wait..." << endl;
  RF.grow_forest(nPerms,pValues,importanceValues);
  //cout << "Time elapsed: " << float(clock() - time_start)/CLOCKS_PER_SEC << " seconds" << endl;
  
  vector<size_t> refIcs(treedata.nfeatures());
  //vector<string> fnames = treedata.featureheaders();
  datadefs::sort_and_make_ref<num_t>(pValues,refIcs);
  datadefs::sort_from_ref<num_t>(importanceValues,refIcs);
  //datadefs::sort_from_ref<string>(fnames,ref_ics);
  
  string targetHeader = treedata.get_featureheader(targetIdx);

  if(pValues[0] <= pValueThreshold)
    {
      //cout << "Writing associations to file '" << output << "'" << endl;
      //ofstream os(output.c_str());
      FILE* po;
      po = fopen(output.c_str(),"w");
      
      for(size_t i = 0; i < treedata.nfeatures(); ++i)
	{
	  if(pValues[i] > pValueThreshold)
	    {
	      break;
	    }

	  fprintf(po,"%s\t%s\t%9.8f\t%9.8f\t%9.8f\n",targetHeader.c_str(),treedata.get_featureheader(refIcs[i]).c_str(),pValues[i],importanceValues[i],treedata.corr(targetIdx,refIcs[i]));

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
