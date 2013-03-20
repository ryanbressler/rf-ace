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
    vector<string> rawTrainData = utils::split(nodeMap["DATA"],',');

    // Set the prediction and data for the node
    if ( isTargetNumerical_ || (!isTargetNumerical_ && forestType_ == forest_t::GBT) ) {
      nodep->setNumTrainPrediction( utils::str2<num_t>(rawTrainPrediction) );
      vector<num_t> numTrainData(rawTrainData.size());
      transform(rawTrainData.begin(),rawTrainData.end(),numTrainData.begin(),utils::str2<num_t>);
      nodep->setNumTrainData(numTrainData);
    } else { 
      nodep->setCatTrainPrediction( rawTrainPrediction );
      vector<cat_t> catTrainData(rawTrainData);
      //transform(rawTrainData.begin(),rawTrainData.end(),catTrainData.begin(),utils::str2<cat_t>);
      nodep->setCatTrainData(catTrainData);
    }
    
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

      num_t splitFitness = utils::str2<num_t>(nodeMap["DI"]);
      
      if ( nodeMap["SPLITTERTYPE"] == "NUMERICAL" ) {

        nodep->setSplitter(splitFitness,nodeMap["SPLITTER"], utils::str2<num_t>(nodeMap["LVALUES"]), lChild, rChild);

      } else if ( nodeMap["SPLITTERTYPE"] == "CATEGORICAL" ){

        unordered_set<string> splitLeftValues = utils::keys(nodeMap["LVALUES"], ':');

        nodep->setSplitter(splitFitness,nodeMap["SPLITTER"], splitLeftValues, lChild, rChild);

      } else if ( nodeMap["SPLITTERTYPE"] == "TEXTUAL" ) {
        nodep->setSplitter(splitFitness,nodeMap["SPLITTER"], utils::str2<uint32_t>(nodeMap["LVALUES"]), lChild, rChild);
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

  //featuresInTree_.clear();

  PredictionFunctionType predictionFunctionType;

  if ( trainData->feature(targetIdx)->isNumerical() ) {
    predictionFunctionType = Node::MEAN;
  } else if ( !trainData->feature(targetIdx)->isNumerical() && forestOptions->forestType == forest_t::GBT ) {
    predictionFunctionType = Node::GAMMA;
  } else {
    predictionFunctionType = Node::MODE;
  }

  nLeaves_ = 1;

  size_t nChildren = 0;

  //Start the recursive node splitting from the root node. This will generate the tree.
  this->recursiveNodeSplit(trainData,
			   targetIdx,
			   forestOptions,
			   random,
			   predictionFunctionType,
			   pmf,
			   bootstrapIcs_,
			   &nLeaves_,
			   nChildren,
			   children_,
			   splitCache_);
  
  children_.resize(nChildren);
  
}

unordered_map<string,num_t> RootNode::getDI() {

  unordered_map<string,num_t> DI;

  for ( size_t nodeIdx = 0; nodeIdx < children_.size(); ++nodeIdx ) {
    unordered_map<string,num_t>::iterator it(DI.find(children_[nodeIdx].getSplitter().name));
    if ( it == DI.end() ) {
      DI[ children_[nodeIdx].getSplitter().name ] = children_[nodeIdx].getSplitter().fitness;
    } else {
      DI[ children_[nodeIdx].getSplitter().name ] += children_[nodeIdx].getSplitter().fitness;
    }
  }

  unordered_map<string,num_t>::iterator it(DI.find(this->getSplitter().name));
  if ( it == DI.end() ) {
    DI[ this->getSplitter().name ]  = this->getSplitter().fitness;
  } else {
    DI[ this->getSplitter().name ] += this->getSplitter().fitness;
  }

  return(DI);

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

const Node::Prediction& RootNode::getPrediction(Treedata* testData, const size_t sampleIdx) {
  return( this->percolate(testData,sampleIdx)->getPrediction() );
}


vector<num_t> RootNode::getChildLeafNumTrainData(Treedata* treeData, const size_t sampleIdx) {

  vector<Node*> leaves = this->percolate(treeData,sampleIdx)->getSubTreeLeaves();

  vector<num_t> allTrainData;

  for ( size_t i = 0; i < leaves.size(); ++i ) {
    vector<num_t> trainData = leaves[i]->getPrediction().numTrainData;
    for ( size_t j = 0; j < trainData.size(); ++j ) {
      allTrainData.push_back(trainData[j]);
    }
  }

  return(allTrainData);
  
}

vector<cat_t> RootNode::getChildLeafCatTrainData(Treedata* treeData, const size_t sampleIdx) {

  vector<Node*> leaves = this->percolate(treeData,sampleIdx)->getSubTreeLeaves();

  vector<cat_t> allTrainData;

  for ( size_t i = 0; i < leaves.size(); ++i ) {
    vector<cat_t> trainData = leaves[i]->getPrediction().catTrainData;
    for ( size_t j = 0; j < trainData.size(); ++j ) {
      allTrainData.push_back(trainData[j]);
    }
  }

  return(allTrainData);

}


