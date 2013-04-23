#ifndef NUMPREDICTTHREAD_HPP
#define NUMPREDICTTHREAD_HPP

#include <vector>

#include "rootnode.hpp"
#include "treedata.hpp"

using namespace std;

class NumPredictThread {
public:
  TreeData* testData_;
  const vector<RootNode*>& rootNodes_;
  forest_t forestType_;
  vector<size_t>& sampleIcs_;
  vector<num_t>* predictions_;
  vector<num_t>* confidence_;
  vector<num_t>& GBTConstants_;
  num_t& GBTShrinkage_;

  NumPredictThread(TreeData* pTestData, 
		   vector<RootNode*>& pRootNodes,
		   forest_t pForestType, 
		   vector<size_t>& pSampleIcs,
		   vector<num_t>* pPredictions, 
		   vector<num_t>* pConfidence,
		   vector<num_t>& pGBTConstants, 
		   num_t& pGBTShrinkage);
  
  ~NumPredictThread();

  void operator()();
};

#endif
