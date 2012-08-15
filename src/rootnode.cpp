#include <string>
#include "rootnode.hpp"
#include "datadefs.hpp"

RootNode::RootNode(Treedata* trainData,
		   options::General_options* parameters,
		   size_t threadIdx): 
  Node(parameters,threadIdx),
  //parameters_(parameters),
  //threadIdx_(threadIdx),
  trainData_(trainData),
  nNodes_(1) {

  trainPredictionCache_.clear();
  
  if ( trainData_ ) {
    trainPredictionCache_.resize(trainData_->nSamples(),datadefs::NUM_NAN);    
  } 
  
}

RootNode::~RootNode() { /* EMPTY DESTRUCTOR */ }

void RootNode::growTree() {

  assert( trainData_ );

  parameters_->validateParameters();

  if ( this->hasChildren() ) {
    this->deleteTree();
    nNodes_ = 1;
  }

  if ( false ) {
    cout << "Growing a tree in thread " << threadIdx_ << " : samplewithReplacement=" << parameters_->sampleWithReplacement 
         << " inBoxFraction=" << parameters_->inBoxFraction << " nMaxLeaves=" << parameters_->nMaxLeaves
         << " nodeSize=" << parameters_->nodeSize << " isRandomSplit=" << parameters_->isRandomSplit << " mTry=" << parameters_->mTry << endl; 
  }

  //Generate the vector for bootstrap indices
  vector<size_t> bootstrapIcs;
  
  size_t targetIdx = trainData_->getFeatureIdx( parameters_->targetStr );
  
  if ( targetIdx == trainData_->end() ) {
    cerr << "Missing target, cannot grow trees!" << endl;
    exit(1);
  }

  //Generate bootstrap indices and oob-indices
  trainData_->bootstrapFromRealSamples(parameters_->randIntGens[threadIdx_], parameters_->sampleWithReplacement, parameters_->inBoxFraction, targetIdx, bootstrapIcs_, oobIcs_);

  //This is to check that the bootstrap sample doesn't contain any missing values (it shouldn't!)
  if ( false ) {
    vector<num_t> targetData = trainData_->getFeatureData(targetIdx,bootstrapIcs_);
    for(size_t i = 0; i < targetData.size(); ++i) {
      assert(!datadefs::isNAN(targetData[i]));
    }

    targetData = trainData_->getFeatureData(targetIdx,oobIcs_);
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

  vector<size_t> featureIcs = utils::range( trainData_->nFeatures() );
  featureIcs.erase( featureIcs.begin() + targetIdx );

  PredictionFunctionType predictionFunctionType;

  if ( trainData_->isFeatureNumerical(targetIdx) ) {
    predictionFunctionType = Node::MEAN;
  } else if ( !trainData_->isFeatureNumerical(targetIdx) && parameters_->modelType == options::GBT ) {
    predictionFunctionType = Node::GAMMA;
  } else {
    predictionFunctionType = Node::MODE;
  }

  size_t nLeaves = 1;

  if ( !parameters_->isRandomSplit ) {
    parameters_->mTry = featureIcs.size();
  }

  //Start the recursive node splitting from the root node. This will generate the tree.
  this->recursiveNodeSplit(trainData_,targetIdx,predictionFunctionType,featureIcs,bootstrapIcs_,featuresInTree_,&nLeaves);
  
  nNodes_ = 2 * nLeaves - 1;

}

size_t RootNode::nNodes() {
  assert( trainData_ );
  return( nNodes_ );
}

vector<size_t> RootNode::getOobIcs() {
  assert( trainData_ );
  return( oobIcs_ );
}

size_t RootNode::nOobSamples() {
  assert( trainData_ );
  return( oobIcs_.size() ); 
}

num_t RootNode::getTestPrediction(Treedata* testData, const size_t sampleIdx) {
  
  return( this->percolate(testData,sampleIdx)->getTrainPrediction() );
  
}

string RootNode::getRawTestPrediction(Treedata* testData, const size_t sampleIdx) {

  return( this->percolate(testData,sampleIdx)->getRawTrainPrediction() );

}

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

num_t RootNode::getPermutedTrainPrediction(const size_t featureIdx, 
					   const size_t sampleIdx) {

  assert( trainData_ );

  // If we have the feature in the tree...
  if ( featuresInTree_.find(featureIdx) != featuresInTree_.end() ) {

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
