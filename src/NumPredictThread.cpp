#include "NumPredictThread.hpp"

#include <thread>

#include "math.hpp"

NumPredictThread::NumPredictThread(TreeData* pTestData, 
				   vector<RootNode*>& pRootNodes,
				   forest_t pForestType, 
				   vector<size_t>& pSampleIcs,
				   vector<num_t>* pPredictions, 
				   vector<num_t>* pConfidence,
				   vector<num_t>& pGBTConstants, 
				   num_t& pGBTShrinkage) :
  testData_(pTestData),
  rootNodes_(pRootNodes),
  forestType_(pForestType),
  sampleIcs_(pSampleIcs),
  predictions_(pPredictions),
  confidence_(pConfidence),
  GBTConstants_(pGBTConstants),
  GBTShrinkage_(pGBTShrinkage)
{ }

NumPredictThread::~NumPredictThread() { }

void NumPredictThread::operator()() {

  size_t nTrees = rootNodes_.size();
  for (size_t i = 0; i < sampleIcs_.size(); ++i) {
    size_t sampleIdx = sampleIcs_.at(i);
    vector<num_t> predictionVec(nTrees);
    for (size_t treeIdx = 0; treeIdx < nTrees; ++treeIdx) {
      predictionVec.at(treeIdx) = rootNodes_.at(treeIdx)->getPrediction(testData_,sampleIdx).numTrainPrediction;
    }
    if (forestType_ == forest_t::GBT) {
      (*predictions_).at(sampleIdx) = GBTConstants_[0];
      (*confidence_).at(sampleIdx) = 0.0;
      for (size_t treeIdx = 0; treeIdx < nTrees; ++treeIdx) {
	predictionVec.at(treeIdx) *= GBTShrinkage_;
        (*predictions_).at(sampleIdx) += predictionVec.at(treeIdx);
      }
    } else {
      (*predictions_).at(sampleIdx) = math::mean(predictionVec);
    }
    (*confidence_).at(sampleIdx)  = sqrt(math::var(predictionVec));
  }
}
