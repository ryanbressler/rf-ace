#include <string>
#include <cmath>
#include <stack>
#include <unordered_set>
#include "math.hpp"
#include "rootnode.hpp"
#include "datadefs.hpp"

RootNode::RootNode():
  children_(0),
  nNodes_(1),
  bootstrapIcs_(0),
  oobIcs_(0),
  minDistToRoot_(0) { /* EMPTY CONSTRUCTOR */ }

RootNode::~RootNode() { /* EMPTY DESTRUCTOR */ }

size_t RootNode::getTreeSizeEstimate(const size_t nSamples, const size_t nMaxLeaves, const size_t nodeSize) const {

  // Upper bound for the number of nodes as dictated by nMaxLeaves, 
  // assuming ternary splits (left,right,missing)
  size_t S1 = static_cast<size_t>( powf(3,ceil(logf(nMaxLeaves-1)/logf(2))) + 1 );

  // Upper bound for the depth of tree as dictated by nSamples, 
  // assuming balanced ternary splits (left,right,missing) 
  size_t k = static_cast<size_t>( ceil( ( logf(nSamples) - logf(nodeSize) ) / logf(3) ) );

  // Tree depth converted to the number of nodes in the tree 
  // S = 3^1 + 3^2 + ... + 3^(k+1)
  size_t S2 = static_cast<size_t>( ( powf(3,k+1) - 3 ) / 2 );

  // Return the smaller of the two upper bounds, S1 and S2
  return( S1 < S2 ? S1 : S2 );

}

void RootNode::reset(const size_t nNodes) {

  assert(nNodes > 1);

  children_.clear();
  children_.resize(nNodes-1);

}

void RootNode::growTree(Treedata* trainData, const size_t targetIdx, const distributions::PMF* pmf, const ForestOptions* forestOptions, distributions::Random* random) {

  size_t nChildren = this->getTreeSizeEstimate(trainData->nSamples(),forestOptions->nMaxLeaves,forestOptions->nodeSize);

  cout << "RootNode::growTree() -- estimating upper bound for tree size: " << nChildren+1 << endl;

  this->reset(nChildren+1);

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

  nChildren = 0;

  //Start the recursive node splitting from the root node. This will generate the tree.
  this->recursiveNodeSplit(trainData,
			   targetIdx,
			   forestOptions,
			   random,
			   predictionFunctionType,
			   pmf,
			   bootstrapIcs_,
			   treeDist,
			   featuresInTree_,
			   minDistToRoot_,
			   &nLeaves,
			   nChildren,
			   children_);
  
  // nNodes_ = 2 * nLeaves - 1;
  nNodes_ = nChildren + 1;

  children_.resize(nChildren);

  this->verifyIntegrity();

}

vector<pair<size_t,size_t> > RootNode::getMinDistFeatures() {

  vector<pair<size_t,size_t> > minDistFeatures;

  for ( set<size_t>::const_iterator it( featuresInTree_.begin() ); it != featuresInTree_.end(); ++it ) {
    size_t featureIdx = *it;
    minDistFeatures.push_back( pair<size_t,size_t>(featureIdx,minDistToRoot_[featureIdx]) );
  }

  return(minDistFeatures);

}

/*
 * - no referring back to root
 * - no double referral to same child
 * - all children referred
 */

void RootNode::verifyIntegrity() const {

  size_t nNodes = this->nNodes();

  assert( children_.size() == nNodes - 1 );

  stack<const Node*> nodesToVisit;
  nodesToVisit.push(this);

  unordered_set<const Node*> nodesReferred;
  nodesReferred.reserve(nNodes);

  while ( ! nodesToVisit.empty() ) {

    const Node* node = nodesToVisit.top();
    nodesToVisit.pop();

    assert( nodesReferred.find(node) == nodesReferred.end() );

    nodesReferred.insert(node);

    if ( node->hasChildren() ) { 
      nodesToVisit.push(node->leftChild());
      nodesToVisit.push(node->rightChild());
    }
    
    if ( node->missingChild() ) {
      nodesToVisit.push(node->missingChild());
    }
    
  }
  
  if ( nodesReferred.size() != nNodes ) {
    cerr << "RootNode::verifyIntegrity() -- only " << nodesReferred.size() << " / " << nNodes << " nodes are reachable!" << endl;
    exit(1);
  }
  
}

// This is a bad function, it exposes the private data to public!!
Node& RootNode::childRef(const size_t childIdx) {

  assert( childIdx < children_.size() );

  return( children_[childIdx] );

}

size_t RootNode::nNodes() const {
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

