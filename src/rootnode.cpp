#include <string>
#include "rootnode.hpp"
#include "datadefs.hpp"

RootNode::RootNode(Treedata* treeData,
		   const size_t targetIdx): 
  Node(),
  treeData_(treeData),
  targetIdx_(targetIdx),
  nNodes_(1),
  trainPredictionCache_(treeData_->nSamples(),datadefs::NUM_NAN) { /* EMPTY CONSTRUCTOR */ }

RootNode::~RootNode() { /* EMPTY DESTRUCTOR */ }

void RootNode::growTree(const GrowInstructions& GI) {

  GI.validate();

  if ( this->hasChildren() ) {
    this->deleteTree();
    nNodes_ = 1;
  }

  if ( false ) {
    cout << "Growing a tree: samplewithReplacement=" << GI.sampleWithReplacement 
         << " sampleSizeFraction=" << GI.sampleSizeFraction << " maxNodesToStop=" << GI.maxNodesToStop
         << " minNodeSizeToStop=" << GI.minNodeSizeToStop << " isRandomSplit=" << GI.isRandomSplit << " nFeaturesForSplit=" << GI.nFeaturesForSplit << endl; 
  }

  //Generate the vector for bootstrap indices
  vector<size_t> bootstrapIcs;
  
  //Generate bootstrap indices and oob-indices
  treeData_->bootstrapFromRealSamples(GI.sampleWithReplacement, GI.sampleSizeFraction, targetIdx_, bootstrapIcs_, oobIcs_);

  //This is to check that the bootstrap sample doesn't contain any missing values (it shouldn't!)
  if(false) {
    vector<num_t> targetData = treeData_->getFeatureData(targetIdx_,bootstrapIcs_);
    for(size_t i = 0; i < targetData.size(); ++i) {
      assert(!datadefs::isNAN(targetData[i]));
    }

    targetData = treeData_->getFeatureData(targetIdx_,oobIcs_);
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

  //Start the recursive node splitting from the root node. This will generate the tree.
  this->recursiveNodeSplit(treeData_,targetIdx_,bootstrapIcs_,GI,featuresInTree_,&nNodes_);
  

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

Node* RootNode::percolateSampleIdx(const size_t sampleIdx) {

  Node* nodep( this );

  // Keep percolating until we hit the leaf
  while ( nodep->hasChildren() ) {

    // Get the splitter feature index
    size_t featureIdxNew = nodep->splitterIdx();

    // Get the respective sample of the splitter feature
    num_t value = treeData_->getFeatureData(featureIdxNew,sampleIdx);

    Node* childNode;

    // Precolate the value, and as a result get a pointer to a child node
    if ( treeData_->isFeatureNumerical(featureIdxNew) ) {
      childNode = nodep->percolateData(value);
    } else {
      childNode = nodep->percolateData(treeData_->getRawFeatureData(featureIdxNew,value));
    }

    // 
    if ( !childNode ) {
      //break;
      num_t r = treeData_->getRandomUnif();
      if ( r <= nodep->leftFraction() ) {
	childNode = nodep->leftChild();
      } else {
	childNode = nodep->rightChild();
      }
    }

    // Update the pointer and continue percolating
    nodep = childNode;
  }

  return( nodep );

}

Node* RootNode::percolateSampleIdx(Treedata* testData, const size_t sampleIdx) {

  Node* nodep( this );

  // Keep percolating until we hit the leaf
  while ( nodep->hasChildren() ) {

    // Get the splitter feature index
    size_t featureIdxNew = testData->getFeatureIdx( nodep->splitterName() );

    Node* childNode = NULL;

    // If the feature exists in the matrix...
    if ( featureIdxNew != testData->end() ) {
      
      // Get the respective sample of the splitter feature
      num_t value = testData->getFeatureData(featureIdxNew,sampleIdx);
      
      // Precolate the value, and as a result get a pointer to a child node
      if ( testData->isFeatureNumerical(featureIdxNew) ) {
	childNode = nodep->percolateData(value);
      } else {
	childNode = nodep->percolateData(testData->getRawFeatureData(featureIdxNew,value));
      }

    }
    
    // The percolation of data to the child nodes failed due to:
    // 1. no feature to split with in the matrix
    // 2. the value to split with is NaN
    // in which case we end up here
    if ( !childNode ) {
      //cout << "Flipping coin..." << endl;
      num_t r = treeData_->getRandomUnif();
      if ( r <= nodep->leftFraction() ) {
        childNode = nodep->leftChild();
      } else {
        childNode = nodep->rightChild();
      }
    }

    // Update the pointer and continue percolating
    nodep = childNode;
  }

  return( nodep );

}

Node* RootNode::percolateSampleIdxAtRandom(const size_t featureIdx, const size_t sampleIdx) {

  Node* nodep( this );

  while ( nodep->hasChildren() ) {

    size_t featureIdxNew = nodep->splitterIdx();

    num_t value = datadefs::NUM_NAN;

    if(featureIdx == featureIdxNew) {
      treeData_->getRandomData(featureIdxNew,value);
    } else {
      value = treeData_->getFeatureData(featureIdxNew,sampleIdx);
    }

    Node* childNode;

    if ( treeData_->isFeatureNumerical(featureIdxNew) ) {
      childNode = nodep->percolateData(value);
    } else {
      childNode = nodep->percolateData(treeData_->getRawFeatureData(featureIdxNew,value));
    }

    if ( !childNode ) {
      //break;
      num_t r = treeData_->getRandomUnif();
      if ( r <= nodep->leftFraction() ) {
        childNode = nodep->leftChild();
      } else {
        childNode = nodep->rightChild();
      }
    }

    nodep = childNode;

  }

  return( nodep );

}


num_t RootNode::getTestPrediction(Treedata* testData, const size_t sampleIdx) {

  return( this->percolateSampleIdx(testData,sampleIdx)->getTrainPrediction() );

}

string RootNode::getRawTestPrediction(Treedata* testData, const size_t sampleIdx) {

  return( this->percolateSampleIdx(testData,sampleIdx)->getRawTrainPrediction() );

}

num_t RootNode::getTrainPrediction(const size_t sampleIdx) {

  // If we don't yet have a prediction for the sample index...
  if ( datadefs::isNAN(trainPredictionCache_[sampleIdx]) ) {

    // We make the prediction!
    trainPredictionCache_[sampleIdx] = this->percolateSampleIdx(sampleIdx)->getTrainPrediction();

    //if ( trainPredictionCache_[sampleIdx] > 1e10 ) {
    //  cout << " " << trainPredictionCache_[sampleIdx];
    //}

  }

  //if ( trainPredictionCache_[sampleIdx] > 1e10 ) {
  //  cout << " " << trainPredictionCache_[sampleIdx];
  //}

  // Return prediction from the cache
  return( trainPredictionCache_[sampleIdx] );

}

num_t RootNode::getPermutedTrainPrediction(const size_t featureIdx, 
					 const size_t sampleIdx) {

  // If we have the feature in the tree...
  if ( featuresInTree_.find(featureIdx) != featuresInTree_.end() ) {

    // We make the prediction by permuting the splitter with the feature!
    return( this->percolateSampleIdxAtRandom(featureIdx,sampleIdx)->getTrainPrediction() );
  
  }
  
  // Otherwise make a regular prediction
  return( this->getTrainPrediction(sampleIdx) );
  
}

vector<num_t> RootNode::getTrainPrediction() {

  size_t nSamples = treeData_->nSamples();

  vector<num_t> prediction(nSamples);

  // predict for all samples
  for ( size_t sampleIdx = 0; sampleIdx < nSamples; ++sampleIdx) {
    prediction[sampleIdx] = this->getTrainPrediction(sampleIdx);
    // cout << "Sample " << i << ", prediction " << curPrediction[i]  << endl;
  }

  return( prediction );
}
