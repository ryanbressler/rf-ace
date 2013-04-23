#include "TreeGrowerThread.hpp"

#include <thread>

TreeGrowerThread::TreeGrowerThread(vector<RootNode*>& pRootNodes,
				   TreeData* pTrainData,
				   const size_t pTargetIdx,
				   const ForestOptions* pForestOptions,
				   const distributions::PMF* pPMF,
				   distributions::Random* pRandomDist) :
  rootNodes_(pRootNodes),
  trainData_(pTrainData),
  targetIdx_(pTargetIdx),
  forestOptions_(pForestOptions),
  PMF_(pPMF),
  randomDist_(pRandomDist)
{ }

TreeGrowerThread::~TreeGrowerThread() {}

void TreeGrowerThread::operator()() {
  for (size_t i = 0; i < rootNodes_.size(); ++i) {
    rootNodes_.at(i)->growTree(trainData_,
				targetIdx_,
				PMF_,
				forestOptions_,
				randomDist_);
  }
}
