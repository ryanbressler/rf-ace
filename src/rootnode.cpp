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


/*
  map<Node*,vector<size_t> > RootNode::percolateSampleIcs(const vector<size_t>& sampleIcs) {
  
  // This horrible construct stores information about which
  // nodes in the tree contain which training samples
  map<Node*,vector<size_t> > percolatedSampleIcs;
  
  // Loop through all train sample indices
  for ( size_t i = 0; i < sampleIcs.size(); ++i) {
  
  Node* nodep( this->percolateSampleIdx(sampleIcs[i]) );
  
  percolatedSampleIcs[nodep].push_back(sampleIcs[i]);
  
  }
  
  if(false) {
  cout << "Train samples percolated accordingly:" << endl;
  size_t iter = 0;
  for(map<Node*,vector<size_t> >::const_iterator it(percolatedSampleIcs.begin()); it != percolatedSampleIcs.end(); ++it, ++iter) {
  cout << "leaf node " << iter << " -- prediction " << it->first->getTrainPrediction() << " :";
  for(size_t i = 0; i < it->second.size(); ++i) {
  cout << " " << it->second[i];
  }
  cout << endl;
  }
  }
  return( percolatedSampleIcs );
  }
*/

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
      break;
      /*
	num_t r = treeData_->getRandomUnif();
	if ( r <= nodep->leftFraction() ) {
	childNode = nodep->leftChild();
	} else {
	childNode = nodep->rightChild();
	}
      */
    }

    // Update the pointer and continue percolating
    nodep = childNode;
  }

  return( nodep );

}

/*
  map<Node*,vector<size_t> > RootNode::percolateSampleIcsAtRandom(const size_t featureIdx, const vector<size_t>& sampleIcs) {
  
  // This horrible construct stores information about which
  // nodes in the tree contain which training samples
  map<Node*,vector<size_t> > percolatedSampleIcs;
  
  // Loop through all train sample indices
  for ( size_t i = 0; i < sampleIcs.size(); ++i) {
  
  // Initialize a pointer to the root node
  Node* nodep( this->percolateSampleIdxAtRandom(featureIdx,sampleIcs[i]) );
  
  percolatedSampleIcs[nodep].push_back(sampleIcs[i]);
  
  }
  
  return( percolatedSampleIcs );
  
  }
*/

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
      break;
    }

    nodep = childNode;

  }

  return( nodep );

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



/*
  num_t RootNode::getPredictionError(const map<Node*,vector<size_t> >& percolatedSampleIcs) {
  
  // Initialize predictionError to 0.0
  num_t predictionError = 0.0;
  
  // Count total number of samples
  //size_t nSamples = 0;
  
  // Get the type of the target
  bool isTargetNumerical = treeData_->isFeatureNumerical(targetIdx_);
  
  // Loop through all nodes that have samples assigned to them
  for ( map<Node*,vector<size_t> >::const_iterator it( percolatedSampleIcs.begin() ); it != percolatedSampleIcs.end(); ++it ) {
  
  // Get target data
  // NOTE1: it->second points to the data indices
  // NOTE2: we happen to know that the indices do not point to data with NANs, which is why
  //        we don't have to make explicit NAN-checks either
  vector<num_t> targetData = treeData_->getFeatureData(targetIdx_,it->second);
  
  // Get node prediction
  // NOTE: it->first points to the node
  num_t nodePrediction = it->first->getTrainPrediction();
  
  assert( !datadefs::isNAN(nodePrediction) );
  
  // Number of percolated sample in node
  size_t nSamplesInNode = targetData.size();
  
  // Depending on the type of the target, different error calculation equation is used
  if ( isTargetNumerical ) {
  
  // Use squared error formula
  for ( size_t i = 0; i < nSamplesInNode; ++i ) {
  predictionError += pow(nodePrediction - targetData[i],2);
  }
  
  } else {
  
  // Use fraction of misprediction formula
  for ( size_t i = 0; i < nSamplesInNode; ++i ) {
  predictionError += nodePrediction != targetData[i];
  }
  }
  
  assert( !datadefs::isNAN(predictionError) );
  
  // Accumulate total sample counter
  //nSamples += nSamplesInNode;
  
  }
  
  // Get mean
  if(nSamples > 0) {
  predictionError /= nSamples;
  } else {
  predictionError = datadefs::NUM_NAN;
  }
  
  
  return( predictionError );
  
  }
*/

/*
  num_t RootNode::getImportance(const size_t featureIdx) {
  
  map<Node*,vector<size_t> > randomlyPercolatedOobSampleIcs = this->percolateSampleIcsAtRandom(featureIdx,oobIcs_);
  
  num_t randomOobPredictionError = this->getPredictionError(randomlyPercolatedOobSampleIcs);
  
  return( randomOobPredictionError - oobError_ );
  
  }
*/
