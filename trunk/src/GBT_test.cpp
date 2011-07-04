
#include<cstdlib>
#include "getopt_pp.h"
#include "GBT.hpp"
//#include "datadefs.hpp"
#include "treedata.hpp"

using namespace std;
using namespace GetOpt;


const size_t DEFAULT_TARGETIDX = 0;
const size_t DEFAULT_NTREES = 0;
//const size_t DEFAULT_MTRY = 0;
const size_t DEFAULT_NODESIZE = 5;
//const size_t DEFAULT_NPERMS = 20;
//const num_t DEFAULT_PTHRESHOLD = 0.10;
//const num_t DEFAULT_ALPHA = 0.95;


int main(int argc, char* argv[])
{

	//------------------------------------------------------------------------
	// 0: parameters
	if(argc == 1 || argc == 2)
	{
		if(argc == 2)
		{
		  string helphandle(argv[1]);
		  if (helphandle != "-h" && helphandle != "--help")
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
		cout << "-n / --ntrees       number of trees per GBT forest (default nsamples)" << endl;
//		cout << "-m / --mtry         number of randomly drawn features per node split (default sqrt(nfeatures))" << endl;
		cout << "-s / --nodesize     minimum number of train samples per node, affects tree depth (default " << DEFAULT_NODESIZE << ")" << endl;
//		cout << "-p / --nperms       number of Random Forests (default " << DEFAULT_NPERMS << ")" << endl;
//		cout << "-t / --pthreshold   p-value threshold below which associations are listed (default " << DEFAULT_PTHRESHOLD << ")" << endl;
//		cout << "-a / --alpha        percentile of contrast importances, defines stringency of the t-test (default " << DEFAULT_ALPHA << ")" << endl;
		cout << endl;
		return EXIT_SUCCESS;
	}

	cout << endl;
	cout << "  ----------------------------------" << endl;
	cout << "  ---  GBT_test version 0.0.2    ---" << endl;
	cout << "  ----------------------------------" << endl;

	//using namespace GetOpt;
	string input = "";
	size_t targetidx = DEFAULT_TARGETIDX;
	size_t ntrees = DEFAULT_NTREES;
//	size_t mtry = DEFAULT_MTRY;
	size_t nodesize = DEFAULT_NODESIZE;
//	size_t nperms = DEFAULT_NPERMS;
//	num_t pthreshold = DEFAULT_PTHRESHOLD;
//	num_t alpha = DEFAULT_ALPHA;
	string output = "";

	GetOpt_pp ops(argc, argv);
		
	ops >> Option('I',"input",input);
	ops >> Option('i',"targetidx",targetidx);
	ops >> Option('n',"ntrees",ntrees);
//	ops >> Option('m',"mtry",mtry);
//	ops >> Option('p',"nperms",nperms);
//	ops >> Option('t',"pthreshold",pthreshold);
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


	//------------------------------------------------------------------------
	// 1: read data into Treedata class (features are rows)
	cout <<endl<< "READ:" << endl;
	Treedata treeData(input);


	//------------------------------------------------------------------------
	// 2: construct a GBT object
	cout <<endl<< "CONSTRUCT:" << endl;
	size_t nTrees(ntrees);
	size_t nMaxLeaves(nodesize);
	float shrinkage(0.5);
	float subSampleSize(0.5);
	GBT gbt(&treeData, nTrees, nMaxLeaves, shrinkage, subSampleSize);

	//------------------------------------------------------------------------
	// 3: select target and grow the forest
	cout <<endl<< "GROWING:" << endl;
	gbt.setTarget(targetidx);
	gbt.allocateForest();
	gbt.growForest();

	//------------------------------------------------------------------------
	// 4: predict using the forest
	cout <<endl<< "PREDICTION:" << endl;
	vector<num_t> prediction(treeData.nSamples());
	gbt.predictForest(&treeData, prediction);

	return(EXIT_SUCCESS);
}
