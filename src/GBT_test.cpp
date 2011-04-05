/*
 * GBT_test.cpp
 *
 *  Created on: Mar 13, 2011
 */

#include<cstdlib>
#include "GBT.hpp"
//#include "datadefs.hpp"
#include "treedata.hpp"

using namespace std;

int main()
{
	// 1: read data into Treedata class (features are rows)
	cout <<endl<< "READ:" << endl;
	bool is_featurerows = true;
	//string fname = "../svn/trunk/data/test_6by10_featurerows_matrix.tsv";
	string fname = "./6_num_features_X_10_cases.tsv";
	Treedata treeData(fname,is_featurerows);


	// 2: construct a GBT object
	cout <<endl<< "CONSTRUCT:" << endl;
	size_t nTrees(5);
	size_t nMaxLeaves(9);
	float shrinkage(0.2);
	float subSampleSize(0.5);
	GBT gbt(&treeData, nTrees, nMaxLeaves, shrinkage, subSampleSize);

	// 3: select target and grow the forest
	cout <<endl<< "GROWING:" << endl;
	gbt.selectTarget(0);
	gbt.growForest();

	// 4: predict using the forest
	cout <<endl<< "PREDICTION:" << endl;
	vector<num_t> prediction(treeData.nsamples());
	gbt.predictForest(&treeData, prediction);

	return(EXIT_SUCCESS);
}
