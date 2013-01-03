#include <string>
#include <cmath>
#include <stack>
#include <unordered_set>
#include "math.hpp"
#include "rootnode.hpp"
#include "datadefs.hpp"

using datadefs::forest_t;

RootNode::RootNode():
  children_(0),
  nLeaves_(0),
  bootstrapIcs_(0),
  oobIcs_(0),
  minDistToRoot_(0) { /* EMPTY CONSTRUCTOR */ }

RootNode::~RootNode() { /* EMPTY DESTRUCTOR */ }

size_t RootNode::getTreeSizeEstimate(const size_t nSamples, const size_t nMaxLeaves, const size_t nodeSize) const {

  // Upper bound for the number of nodes as dictated by nMaxLeaves, 
  // assuming ternary splits (left,right,missing)
  size_t S1 = static_cast<size_t>(ceil(1.0 + 3.0 * ( nMaxLeaves - 1 ) / 2.0)); //static_cast<size_t>( powf(3,ceil(logf(nMaxLeaves-1)/logf(2))) + 1 );

  // Upper bound for the number of nodes as dictated by sample size,
  // assuming ternary splits (left,right,missing)
  size_t S2 = static_cast<size_t>(ceil(1.0 + 3.0 * ( ceil(nSamples/nodeSize) - 1 ) / 2.0 ));

  // Upper bound for the depth of tree as dictated by nSamples, 
  // assuming balanced ternary splits (left,right,missing) 
  // size_t k = static_cast<size_t>( ceil( ( logf(nSamples) - logf(nodeSize) ) / logf(3) ) );

  // Tree depth converted to the number of nodes in the tree 
  // S = 1 + 3^1 + 3^2 + ... + 3^(k+1)
  // size_t S2 = 1 + static_cast<size_t>( ( powf(3,k+1) - 3 ) / 2 );

  //cout << "f(" << nSamples << "," << nMaxLeaves << "," << nodeSize << ") = min(" << S1 << "," << S2 << ")" << endl;  

  // Return the smaller of the two upper bounds, S1 and S2
  return( S1 < S2 ? S1 : S2 );

}

void RootNode::reset(const size_t nNodes) {

  assert(nNodes > 1);

  children_.clear();
  children_.resize(nNodes-1);

}

void RootNode::growTree(Treedata* trainData, const size_t targetIdx, const distributions::PMF* pmf, const ForestOptions* forestOptions, distributions::Random* random) {

  size_t nMaxNodes = this->getTreeSizeEstimate(trainData->nRealSamples(targetIdx),forestOptions->nMaxLeaves,forestOptions->nodeSize);

  this->reset(nMaxNodes);

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
  } else if ( !trainData->isFeatureNumerical(targetIdx) && forestOptions->forestType == forest_t::GBT ) {
    predictionFunctionType = Node::GAMMA;
  } else {
    predictionFunctionType = Node::MODE;
  }

  nLeaves_ = 1;

  size_t treeDist = 0;

  minDistToRoot_.clear();
  minDistToRoot_.resize(2*trainData->nFeatures(),datadefs::MAX_IDX);

  size_t nChildren = 0;

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
			   &nLeaves_,
			   nChildren,
			   children_);
  
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

  stack<const Node*> nodesToVisit;
  nodesToVisit.push(this);

  unordered_set<const Node*> nodesReferred;
  nodesReferred.rehash(4*nNodes);

  while ( ! nodesToVisit.empty() ) {

    const Node* node = nodesToVisit.top();
    nodesToVisit.pop();

    if ( nodesReferred.find(node) != nodesReferred.end() ) {
      cerr << "RootNode::verifyIntegrity() -- double referral to the same node in the tree, which should be impossible!" << endl;
      exit(1);
    }

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
  return( children_.size() + 1 );
}

size_t RootNode::nLeaves() const {

  if ( nLeaves_ == 0 ) {
    cerr << "ERROR: RootNode::nLeaves() -- trying to save a forest that has been read from file? Don't!" << endl;
    exit(1);
  }

  return( nLeaves_ );
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

