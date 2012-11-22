#include <string>
#include "math.hpp"
#include "rootnode.hpp"
#include "datadefs.hpp"

RootNode::RootNode(): 
  nNodes_(1),
  bootstrapIcs_(0),
  oobIcs_(0),
  minDistToRoot_(0) { /* EMPTY CONSTRUCTOR */ }

RootNode::~RootNode() { /* EMPTY DESTRUCTOR */ }

void RootNode::growTree(Treedata* trainData, const size_t targetIdx, const distributions::PMF* pmf, const ForestOptions* forestOptions, distributions::Random* random) {

  if ( this->hasChildren() ) {
    this->deleteTree();
    nNodes_ = 1;
  }

  //Generate the vector for bootstrap indices
  vector<size_t> bootstrapIcs;
  
  //Generate bootstrap indices and oob-indices
  trainData->bootstrapFromRealSamples(random, forestOptions->sampleWithReplacement, forestOptions->inBoxFraction, targetIdx, bootstrapIcs_, oobIcs_);

  //This is to check that the bootstrap sample doesn't contain any missing values (it shouldn't!)
  if ( false ) {
    vector<num_t> targetData = trainData->getFeatureData(targetIdx,bootstrapIcs_);
    for(size_t i = 0; i < targetData.size(); ++i) {
      assert(!datadefs::isNAN(targetData[i]));
    }

    targetData = trainData->getFeatureData(targetIdx,oobIcs_);
    for(size_t i = 0; i < targetData.size(); ++i) {
      assert(!datadefs::isNAN(targetData[i]));
    }
    cout << "bootstrap samples look ok, no missing values detected" << endl;
  }

  if(false) {
    cout << "tree bootstrap indices [";
    for(size_t i = 0; i < bootstrapIcs_.size(); ++i) {
      cout << " " << bootstrapIcs_[i];
    }
    cout << " ]  oob [";
    for(size_t i = 0; i < oobIcs_.size(); ++i) {
      cout << " " << oobIcs_[i];
    }
    cout << " ]" << endl << endl;
  }

  featuresInTree_.clear();

  PredictionFunctionType predictionFunctionType;

  if ( trainData->isFeatureNumerical(targetIdx) ) {
    predictionFunctionType = Node::MEAN;
  } else if ( !trainData->isFeatureNumerical(targetIdx) && forestOptions->forestType == ForestOptions::ForestType::GBT ) {
    predictionFunctionType = Node::GAMMA;
  } else {
    predictionFunctionType = Node::MODE;
  }

  size_t nLeaves = 1;

  size_t treeDist = 0;

  minDistToRoot_.clear();
  minDistToRoot_.resize(2*trainData->nFeatures(),datadefs::MAX_IDX);

  //Start the recursive node splitting from the root node. This will generate the tree.
  this->recursiveNodeSplit(trainData,targetIdx,forestOptions,random,predictionFunctionType,pmf,bootstrapIcs_,treeDist,featuresInTree_,minDistToRoot_,&nLeaves);
  
  nNodes_ = 2 * nLeaves - 1;

}

vector<pair<size_t,size_t> > RootNode::getMinDistFeatures() {

  vector<pair<size_t,size_t> > minDistFeatures;

  for ( set<size_t>::const_iterator it( featuresInTree_.begin() ); it != featuresInTree_.end(); ++it ) {
    size_t featureIdx = *it;
    minDistFeatures.push_back( pair<size_t,size_t>(featureIdx,minDistToRoot_[featureIdx]) );
  }

  return(minDistFeatures);

}

size_t RootNode::nNodes() {
  return( nNodes_ );
}

vector<size_t> RootNode::getOobIcs() {
  return( oobIcs_ );
}

size_t RootNode::nOobSamples() {
  return( oobIcs_.size() ); 
}

num_t RootNode::getTestPrediction(Treedata* testData, const size_t sampleIdx) {
  
  return( this->percolate(testData,sampleIdx)->getTrainPrediction() );
  
}

string RootNode::getRawTestPrediction(Treedata* testData, const size_t sampleIdx) {

  return( this->percolate(testData,sampleIdx)->getRawTrainPrediction() );

}

/*
  num_t RootNode::getTrainPrediction(const size_t sampleIdx) {
  
  assert( trainData_ );
  
  // If we don't yet have a prediction for the sample index...
  if ( datadefs::isNAN(trainPredictionCache_[sampleIdx]) ) {
  
  // We make the prediction!
  trainPredictionCache_[sampleIdx] = this->percolate(trainData_,sampleIdx)->getTrainPrediction();
  
  }
  
  // Return prediction from the cache
  return( trainPredictionCache_[sampleIdx] );
  
  }
*/

/*
  num_t RootNode::getPermutedTrainPrediction(const size_t featureIdx, 
  const size_t sampleIdx) {
  
  assert( trainData_ );
  
  // If we have the feature in the tree...
  if ( minDistToRoot_.find(featureIdx) != minDistToRoot_.end() ) {
  
  // We make the prediction by permuting the splitter with the feature!
  return( this->percolate(trainData_,sampleIdx,featureIdx)->getTrainPrediction() );
  
  }
  
  // Otherwise make a regular prediction
  return( this->getTrainPrediction(sampleIdx) );
  
  }
  
  vector<num_t> RootNode::getTrainPrediction() {
  
  assert( trainData_ );
  
  size_t nSamples = trainData_->nSamples();
  
  vector<num_t> prediction(nSamples);
  
  // predict for all samples
  for ( size_t sampleIdx = 0; sampleIdx < nSamples; ++sampleIdx) {
  prediction[sampleIdx] = this->getTrainPrediction(sampleIdx);
  // cout << "Sample " << i << ", prediction " << curPrediction[i]  << endl;
  }
  
  return( prediction );
  }
*/

