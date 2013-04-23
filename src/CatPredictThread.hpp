#ifndef CATPREDICTTHREAD_HPP
#define CATPREDICTTHREAD_HPP

#include <vector>

#include "rootnode.hpp"
#include "treedata.hpp"

using namespace std;

class CatPredictThread {
public:
  TreeData* testData_;
  vector<RootNode*>& rootNodes_;
  forest_t forestType_;
  const vector<size_t>& sampleIcs_;
  vector<cat_t>* predictions_;
  vector<num_t>* confidence_;
  vector<cat_t>& categories_;
  vector<num_t>& GBTConstants_;
  num_t& GBTShrinkage_;
  
  CatPredictThread(TreeData* pTestData, 
		   vector<RootNode*>& pRootNodes,
		   forest_t pForestType,
		   vector<size_t>& pSampleIcs, 
		   vector<cat_t>* pPredictions,
		   vector<num_t>* pConfidence, 
		   vector<cat_t>& pCategories,
		   vector<num_t>& pGBTConstants, 
		   num_t& pGBTShrinkage);
  
  ~CatPredictThread();

  void operator()();
};

#endif
