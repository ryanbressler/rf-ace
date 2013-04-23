#include "CatPredictThread.hpp"

#include <thread>

#include "math.hpp"

CatPredictThread::CatPredictThread(TreeData* pTestData, 
				   vector<RootNode*>& pRootNodes,
				   forest_t pForestType,
				   vector<size_t>& pSampleIcs, 
				   vector<cat_t>* pPredictions,
				   vector<num_t>* pConfidence, 
				   vector<cat_t>& pCategories,
				   vector<num_t>& pGBTConstants, 
				   num_t& pGBTShrinkage) :
  testData_(pTestData),
  rootNodes_(pRootNodes),
  forestType_(pForestType),
  sampleIcs_(pSampleIcs),
  predictions_(pPredictions),
  confidence_(pConfidence),
  categories_(pCategories),
  GBTConstants_(pGBTConstants),
  GBTShrinkage_(pGBTShrinkage)
{ }

CatPredictThread::~CatPredictThread() { }

void CatPredictThread::operator()() {
  size_t nTrees = rootNodes_.size();
  for ( size_t i = 0; i < sampleIcs_.size(); ++i ) {

    size_t sampleIdx = sampleIcs_[i];

    if ( forestType_ == forest_t::GBT ) {

      cerr << "Prediction with GBT is turned OFF" << endl;
      exit(1);

      size_t nCategories = categories_.size();
      num_t maxProb = 0;
      size_t maxProbCat = 0;

      for ( size_t categoryIdx = 0; categoryIdx < nCategories; ++categoryIdx ) {

        num_t cumProb = GBTConstants_[categoryIdx];

        for ( size_t iterIdx = 0; iterIdx < nTrees / nCategories; ++iterIdx ) {
          size_t treeIdx = iterIdx * nCategories + categoryIdx;
          cumProb += GBTShrinkage_ * rootNodes_[treeIdx]->getPrediction(testData_, sampleIdx).numTrainPrediction;
        }

        if (cumProb > maxProb) {
          maxProb = cumProb;
          maxProbCat = categoryIdx;
        }
      }

      (*predictions_).at(sampleIdx) = categories_[maxProbCat];
      (*confidence_).at(sampleIdx) = 0.0;

      // MISSING: confidence
    }
    else {

      vector<string> predictionVec(nTrees);
      for ( size_t treeIdx = 0; treeIdx < nTrees; ++treeIdx ) {
        predictionVec.at(treeIdx) = rootNodes_.at(treeIdx)->getPrediction(testData_, sampleIdx).catTrainPrediction;
      }

      (*predictions_).at(sampleIdx) = math::mode(predictionVec);
      (*confidence_).at(sampleIdx) = 1.0 * math::nMismatches(predictionVec, (*predictions_)[sampleIdx]) / nTrees;
    }
  }
}
