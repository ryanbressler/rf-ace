#include <string>
#include <cmath>
#include <stack>
#include <unordered_set>
#include "math.hpp"
#include "rootnode.hpp"
#include "datadefs.hpp"

using datadefs::forest_t;

RootNode::RootNode() {}

RootNode::RootNode(Treedata* trainData, const size_t targetIdx, const distributions::PMF* pmf, const ForestOptions* forestOptions, distributions::Random* random):
  forestType_(forestOptions->forestType),
  targetName_(trainData->feature(targetIdx)->name()),
  isTargetNumerical_(trainData->feature(targetIdx)->isNumerical()),
  children_(0),
  nLeaves_(0),
  bootstrapIcs_(0),
  oobIcs_(0),
  minDistToRoot_(0) {

  this->growTree(trainData,targetIdx,pmf,forestOptions,random);

}

RootNode::RootNode(ifstream& treeStream) {

  this->loadTree(treeStream);

}

RootNode::~RootNode() { /* EMPTY DESTRUCTOR */ }

size_t RootNode::getTreeSizeEstimate(const size_t nSamples, const size_t nMaxLeaves, const size_t nodeSize) const {

  // Upper bound for the number of nodes as dictated by nMaxLeaves, 
  // assuming ternary splits (left,right,missing)
  size_t S1 = static_cast<size_t>(ceil(1.0 + 3.0 * ( nMaxLeaves - 1 ) / 2.0)); 

  // Upper bound for the number of nodes as dictated by sample size,
  // assuming ternary splits (left,right,missing)
  size_t S2 = static_cast<size_t>(ceil(1.0 + 3.0 * ( ceil(nSamples/nodeSize) - 1 ) / 2.0 ));

  // Return the smaller of the two upper bounds, S1 and S2
  return( S1 < S2 ? S1 : S2 );

}

void RootNode::reset(const size_t nNodes) {

  assert(nNodes > 1);

  children_.clear();
  children_.resize(nNodes-1);

}

void RootNode::loadTree(ifstream& treeStream) {

  unordered_map<string,Node*> treeMap;

  size_t nNodes = 0;
  size_t nNodesAllocated = 0;

  string newLine;

  getline(treeStream,newLine);

  // remove trailing end-of-line characters
  newLine = utils::chomp(newLine);

  assert(newLine.compare(0, 5, "TREE=") == 0);
  map<string,string> treeSetup = utils::parse(newLine, ',', '=', '"');
  assert( treeSetup.find("NNODES") != treeSetup.end() );
  nNodes = utils::str2<size_t>(treeSetup["NNODES"]);

  forestType_ = datadefs::forestTypeAssign.at(treeSetup["FOREST"]);
  targetName_ = treeSetup["TARGET"];
  isTargetNumerical_ = utils::str2<bool>(treeSetup["ISTARGETNUMERICAL"]);

  this->reset(nNodes);
  treeMap["*"] = this;

  for ( size_t nodeIdx = 0; nodeIdx < nNodes; ++nodeIdx ) {

    assert( getline(treeStream,newLine) );
    
    // remove trailing end-of-line characters
    newLine = utils::chomp(newLine);

    map<string,string> nodeMap = utils::parse(newLine, ',', '=', '"');
    
    Node* nodep = treeMap[nodeMap["NODE"]];
    
    string rawTrainPrediction = nodeMap["PRED"];
    num_t trainPrediction = datadefs::NUM_NAN;
    
    if ( isTargetNumerical_ || (!isTargetNumerical_ && forestType_ == forest_t::GBT) ) {
      trainPrediction = utils::str2<num_t>(rawTrainPrediction);
    }
    
    vector<string> foo2 = utils::split(nodeMap["DATA"],',');
    vector<num_t> trainData(foo2.size());
    transform(foo2.begin(),foo2.end(),trainData.begin(),utils::str2<num_t>);
    
    // Set train data prediction of the node
    nodep->setTrainPrediction(trainPrediction, rawTrainPrediction);
    nodep->setTrainData(trainData);
    
    // If the node has a splitter... 
    if ( nodeMap.find("SPLITTER") != nodeMap.end() ) {
      
      size_t leftChildIdx = nNodesAllocated++;
      size_t rightChildIdx = nNodesAllocated++;
      
      if ( nNodesAllocated + 1 > nNodes ) {
	cerr << "RootNode::loadForest() -- the tree contains more nodes than declared!" << endl;
	exit(1);
      }
      
      Node& lChild = this->childRef(leftChildIdx);
      Node& rChild = this->childRef(rightChildIdx);
      
      if ( nodeMap["SPLITTERTYPE"] == "NUMERICAL" ) {

        nodep->setSplitter(nodeMap["SPLITTER"], utils::str2<num_t>(nodeMap["LVALUES"]), lChild, rChild);

      } else if ( nodeMap["SPLITTERTYPE"] == "CATEGORICAL" ){

        unordered_set<string> splitLeftValues = utils::keys(nodeMap["LVALUES"], ':');

        nodep->setSplitter(nodeMap["SPLITTER"], splitLeftValues, lChild, rChild);

      } else if ( nodeMap["SPLITTERTYPE"] == "TEXTUAL" ) {
        nodep->setSplitter(nodeMap["SPLITTER"], utils::str2<uint32_t>(nodeMap["LVALUES"]), lChild, rChild);
      } else {
        cerr << "ERROR: incompatible splitter type '" << nodeMap["SPLITTERTYPE"] << endl;
        exit(1);
      }

      assert( &lChild == nodep->leftChild() );
      assert( &rChild == nodep->rightChild() );

      treeMap[nodeMap["NODE"] + "L"] = nodep->leftChild();
      treeMap[nodeMap["NODE"] + "R"] = nodep->rightChild();

      // In case there is a branch for case when the splitter has missing value
      if ( nodeMap.find("M") != nodeMap.end() ) {
        size_t missingChildIdx = nNodesAllocated++;
        Node& mChild = this->childRef(missingChildIdx);
        nodep->setMissingChild(mChild);
        assert( &mChild == nodep->missingChild() );
        treeMap[nodeMap["NODE"] + "M"] = nodep->missingChild();
      }

    }

  }

  assert( nNodesAllocated + 1 == nNodes );

  treeStream.peek();

}

void RootNode::writeTree(ofstream& toFile) {

  toFile << "TREE=," << flush;
  
  if (forestType_ == forest_t::GBT) {
    toFile << "FOREST=GBT";
  } else if (forestType_ == forest_t::RF) {
    toFile << "FOREST=RF";
  } else if (forestType_ == forest_t::QRF) {
    toFile << "FOREST=QRF";
  } else {
    cerr << "StochasticForest::saveForest() -- Unknown forest type!" << endl;
    exit(1);
  }

  toFile << ",NNODES=" << this->nNodes() << ",NLEAVES=" << this->nLeaves() << ",TARGET=\"" << targetName_ << "\",ISTARGETNUMERICAL=" << isTargetNumerical_ << ",CLASS=\"\"" << endl;

  string traversal("*");
  this->recursiveWriteTree(traversal,toFile);

}


void RootNode::growTree(Treedata* trainData, const size_t targetIdx, const distributions::PMF* pmf, const ForestOptions* forestOptions, distributions::Random* random) {

  forestType_ = forestOptions->forestType;
  targetName_ = trainData->feature(targetIdx)->name();
  isTargetNumerical_ = trainData->feature(targetIdx)->isNumerical();

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

  if ( trainData->feature(targetIdx)->isNumerical() ) {
    predictionFunctionType = Node::MEAN;
  } else if ( !trainData->feature(targetIdx)->isNumerical() && forestOptions->forestType == forest_t::GBT ) {
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
			   children_,
			   splitCache_);
  
  children_.resize(nChildren);

  if ( forestOptions->forestType == forest_t::QRF ) {
    //assert( trainData->feature(targetIdx)->isNumerical() );
    for ( size_t i = 0; i < oobIcs_.size(); ++i ) {
      num_t x = trainData->feature(targetIdx)->data[oobIcs_[i]];
      this->percolate(trainData,oobIcs_[i])->addTrainData(x);
    }
  }
  
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

vector<num_t> RootNode::getChildLeafTrainData(Treedata* treeData, const size_t sampleIdx) {

  vector<Node*> leaves = this->percolate(treeData,sampleIdx)->getChildLeaves();

  vector<num_t> allTrainData;

  for ( size_t i = 0; i < leaves.size(); ++i ) {
    vector<num_t> trainData = leaves[i]->getTrainData();
    for ( size_t j = 0; j < trainData.size(); ++j ) {
      allTrainData.push_back(trainData[j]);
    }
  }

  return(allTrainData);
  
}

