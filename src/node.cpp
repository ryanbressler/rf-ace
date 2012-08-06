#include<iostream>
#include<cassert>
#include<iomanip>

#include "node.hpp"
#include "datadefs.hpp"
#include "utils.hpp"
#include "math.hpp"

using namespace std;
using datadefs::num_t;

Node::Node():
  trainPrediction_(datadefs::NUM_NAN),
  rawTrainPrediction_(datadefs::STR_NAN),
  leftChild_(NULL),
  rightChild_(NULL) {
  
}

/** 
 * Destructor for the node: initiates a recursive destruction for all child nodes
 */
Node::~Node() {
  
  //If the node has children, some moery cleanup needs to be performed
  if ( this->hasChildren() ) {

    // Deallocates dynamically allocated memory
    this->deleteTree();
  }

}

/**
 * Deletes child nodes, which will cascade all the way to the leaf nodes 
 */
void Node::deleteTree() {

  if ( false ) {
    cout << "DEL " << leftChild_ << endl;
    cout << "DEL " << rightChild_ << endl; 
  }

  delete leftChild_;
  delete rightChild_;
  
}

// !! Documentation: consider combining with the documentation in the header
// !! file, fleshing it out a bit. Ideally, implementation notes should fall
// !! here; notes on the abstraction should fall in the header file.
void Node::setSplitter(const size_t splitterIdx, 
		       const string& splitterName, 
		       const num_t leftFraction, 
		       const num_t splitLeftLeqValue) {
  
  if ( this->hasChildren() ) {
    cerr << "Cannot set a splitter to a node twice!" << endl;
    exit(1);
  }

  if ( leftFraction < 0.0 || leftFraction > 1.0 ) {
    cerr << "Node::setSplitter() -- leftFraction must be within 0..1!" << endl;
    exit(1);
  }

  splitter_.idx = splitterIdx;
  splitter_.name = splitterName;
  splitter_.isNumerical = true;
  splitter_.leftFraction = leftFraction;
  splitter_.leftLeqValue = splitLeftLeqValue;
  splitter_.leftFraction = leftFraction;

  leftChild_ = new Node();
  rightChild_ = new Node();

  if ( false ) {
    cout << "NEW " << leftChild_ << endl;
    cout << "NEW " << rightChild_ << endl;
  }

}

// !! Documentation: consider combining with the documentation in the header
// !! file, fleshing it out a bit. Ideally, implementation notes should fall
// !! here; notes on the abstraction should fall in the header file.
void Node::setSplitter(const size_t splitterIdx, 
		       const string& splitterName, 
		       const num_t leftFraction, 
		       const set<string>& leftSplitValues, 
		       const set<string>& rightSplitValues) {

  if ( this->hasChildren() ) {
    cerr << "Node::setSplitter() -- cannot set a splitter to a node twice!" << endl;
    exit(1);
  }

  if ( leftFraction < 0.0 || leftFraction > 1.0 ) {
    cerr << "Node::setSplitter() -- leftFraction must be within 0..1!" << endl;
    exit(1);
  }

  splitter_.idx = splitterIdx;
  splitter_.name = splitterName;
  splitter_.isNumerical = false;
  splitter_.leftFraction = leftFraction;
  splitter_.leftValues = leftSplitValues;
  splitter_.rightValues = rightSplitValues;

  leftChild_ = new Node();
  rightChild_ = new Node();

  if ( false ) {
    cout << "NEW " << leftChild_ << endl;
    cout << "NEW " << rightChild_ << endl;
  }

}


Node* Node::percolateData(const num_t data) {
  
  assert( splitter_.isNumerical );

  // Return NULL if the node doesn't have children ( == is a leaf node )
  // or the data is NAN
  if ( !this->hasChildren() || datadefs::isNAN(data) ) {
    return( NULL );
  }

  return( data <= splitter_.leftLeqValue ? leftChild_ : rightChild_ );
   
}

Node* Node::percolateData(const string& data) {

  assert( !splitter_.isNumerical );

  // Return NULL if the node doesn't have children ( == is a leaf node )
  // or the data is NAN
  if ( !this->hasChildren() || datadefs::isNAN_STR(data) ) {
    return( NULL );
  }

  // Return left child if splits left
  if ( splitter_.leftValues.find(data) != splitter_.leftValues.end() ) {
    return( leftChild_ );
  }

  // Return right child if splits right
  if ( splitter_.rightValues.find(data) != splitter_.rightValues.end() ) {
    return( rightChild_ );
  }

  //cout << "Shit, categorical split failed!" << endl;

  // Return this if splits neither left nor right, which can happen if
  // the splitter is categorical
  return( NULL );

}

Node* Node::leftChild() {
  return( leftChild_ );
}


Node* Node::rightChild() {
  return( rightChild_ );
}

/**
 * Recursively prints a tree to a stream (file)
 */
void Node::print(string& traversal, ofstream& toFile) {

  toFile << "NODE=" << traversal << ",PRED=" << rawTrainPrediction_;
  
  if ( this->hasChildren() ) {
    
    string splitterType = splitter_.isNumerical ? "NUMERICAL" : "CATEGORICAL";

    toFile << ",SPLITTER=" << "\"" << splitter_.name << "\""
	   << ",SPLITTERTYPE=" << splitterType
	   << ",LFRACTION=" << splitter_.leftFraction; 

    if ( splitter_.isNumerical ) {
      toFile << ",LVALUES=" << splitter_.leftLeqValue
	     << ",RVALUES=" << splitter_.leftLeqValue << endl;
    } else {

      toFile << ",LVALUES=" << "\""; utils::write(toFile,splitter_.leftValues.begin(),splitter_.leftValues.end(),':'); toFile << "\""
	     << ",RVALUES=" << "\""; utils::write(toFile,splitter_.rightValues.begin(),splitter_.rightValues.end(),':'); toFile << "\"" << endl;
    }
    
    string traversalLeft = traversal;
    traversalLeft.append("L");
    string traversalRight = traversal;
    traversalRight.append("R");
    
    leftChild_->print(traversalLeft,toFile);
    rightChild_->print(traversalRight,toFile);
    
  } else {
    toFile << endl;
  }
}

void Node::print(ofstream& toFile) {
  
  string traversal("*");
  this->print(traversal,toFile);

}


void Node::setTrainPrediction(const num_t trainPrediction, const string& rawTrainPrediction) {
  trainPrediction_ = trainPrediction;
  rawTrainPrediction_ = rawTrainPrediction;
}

// !! Documentation: just your usual accessor, returning a copy of
// !! trainPrediction_.
num_t Node::getTrainPrediction() {
  assert( !datadefs::isNAN(trainPrediction_) );
  return( trainPrediction_ );
}

string Node::getRawTrainPrediction() {
  assert( !datadefs::isNAN_STR(rawTrainPrediction_) );
  return( rawTrainPrediction_ );
}

void Node::recursiveNodeSplit(Treedata* treeData,
			      const size_t targetIdx,
			      options::General_options* parameters,
			      const size_t threadIdx,
			      const PredictionFunctionType& predictionFunctionType,
			      vector<size_t> featureIcs,
			      const vector<size_t>& sampleIcs,
			      set<size_t>& featuresInTree,
			      size_t* nLeaves) {

  if ( false ) {
    cout << "REC " << this << endl;
    cout << " " << treeData << endl;
    cout << " " << parameters << endl;
    cout << " " << nLeaves << endl;
  }

  size_t nSamples = sampleIcs.size();

  assert( *nLeaves <= parameters->nMaxLeaves );

  if ( nSamples < 2 * parameters->nodeSize || *nLeaves == parameters->nMaxLeaves ) {
    
    vector<num_t> trainData = treeData->getFeatureData(targetIdx,sampleIcs);

    if ( predictionFunctionType == MEAN ) {
      num_t trainPrediction = math::mean(trainData);
      string rawTrainPrediction = utils::num2str(trainPrediction);
      this->setTrainPrediction( trainPrediction, rawTrainPrediction );
    } else if ( predictionFunctionType == MODE ) {
      num_t trainPrediction = math::mode<num_t>(trainData);
      string rawTrainPrediction = treeData->getRawFeatureData(targetIdx,trainPrediction);
      this->setTrainPrediction( trainPrediction, rawTrainPrediction );
    } else if ( predictionFunctionType == GAMMA ) {
      num_t trainPrediction = math::gamma(trainData, treeData->nCategories(targetIdx) );
      string rawTrainPrediction = utils::num2str(trainPrediction);
      this->setTrainPrediction( trainPrediction, rawTrainPrediction );
    } else {
      cerr << "Node::recursiveNodeSplit() -- unknown prediction function!" << endl;
      exit(1);
    }

    //assert( !datadefs::isNAN(trainPrediction_) );
    assert( !datadefs::isNAN_STR(rawTrainPrediction_) );

    return;
  }

  vector<size_t> featureSampleIcs = featureIcs;

  if(parameters->isRandomSplit) {

    utils::permute(featureSampleIcs,parameters->randIntGens[threadIdx]);

    // Take only the first ones
    featureSampleIcs.resize(parameters->mTry);

  }

  // With 1% sampling rate assign contrasts
  if(parameters->useContrasts) {
    for(size_t i = 0; i < parameters->mTry; ++i) {
      
      // If the sampled feature is a contrast... 
      if( parameters->randIntGens[threadIdx].uniform() < 0.01 ) { // %1 sampling rate
	
	featureSampleIcs[i] += treeData->nFeatures();
      }
    }
  } 

  //datadefs::print<size_t>(featureSampleIcs);
  
  assert( featureSampleIcs.size() == parameters->mTry );
  
  vector<size_t> sampleIcs_left,sampleIcs_right;

  size_t splitFeatureIdx;
  num_t splitFitness;
  
  bool foundSplit = this->regularSplitterSeek(treeData,
					      targetIdx,
					      sampleIcs,
					      featureSampleIcs,
					      parameters,
					      splitFeatureIdx,
					      sampleIcs_left,
					      sampleIcs_right,
					      splitFitness);
        
  if ( !foundSplit ) {

    vector<num_t> trainData = treeData->getFeatureData(targetIdx,sampleIcs);
    
    if ( predictionFunctionType == MEAN ) {
      num_t trainPrediction = math::mean(trainData);
      string rawTrainPrediction = utils::num2str(trainPrediction);
      this->setTrainPrediction( trainPrediction, rawTrainPrediction );
    } else if ( predictionFunctionType == MODE ) {
      num_t trainPrediction = math::mode<num_t>(trainData);
      string rawTrainPrediction = treeData->getRawFeatureData(targetIdx,trainPrediction);
      this->setTrainPrediction( trainPrediction, rawTrainPrediction );
    } else if ( predictionFunctionType == GAMMA ) {
      num_t trainPrediction = math::gamma(trainData,treeData->nCategories(targetIdx));
      string rawTrainPrediction = utils::num2str(trainPrediction);
      this->setTrainPrediction( trainPrediction, rawTrainPrediction );
    } else {
      cerr << "Node::recursiveNodeSplit() -- unknown prediction function!" << endl;
      exit(1);
    }

    //assert( !datadefs::isNAN(trainPrediction_) );
    assert( !datadefs::isNAN_STR(rawTrainPrediction_) );

    return;
  }
  
  featuresInTree.insert(splitFeatureIdx);
  *nLeaves += 1;
  
  leftChild_->recursiveNodeSplit(treeData,targetIdx,parameters,threadIdx,predictionFunctionType,featureIcs,sampleIcs_left,featuresInTree,nLeaves);
  rightChild_->recursiveNodeSplit(treeData,targetIdx,parameters,threadIdx,predictionFunctionType,featureIcs,sampleIcs_right,featuresInTree,nLeaves);
  
}

bool Node::regularSplitterSeek(Treedata* treeData,
			       const size_t targetIdx,
			       const vector<size_t>& sampleIcs,
			       const vector<size_t>& featureSampleIcs,
			       options::General_options* parameters,
			       size_t& splitFeatureIdx,
			       vector<size_t>& sampleIcs_left,
			       vector<size_t>& sampleIcs_right,
			       num_t& splitFitness) {
  
  // This many features will be tested for splitting the data
  size_t nFeaturesForSplit = featureSampleIcs.size();
  
  num_t splitValue = datadefs::NUM_NAN;
  set<num_t> splitValues_left; 
  set<num_t> splitValues_right;

  // Initialize split fitness to lowest possible value
  splitFitness = 0.0;

  // Loop through candidate splitters
  for ( size_t i = 0; i < nFeaturesForSplit; ++i ) {
    
    // Get split feature index
    size_t newSplitFeatureIdx = featureSampleIcs[i];

    // Get type of the splitter
    bool isFeatureNumerical = treeData->isFeatureNumerical(newSplitFeatureIdx);

    // We don't want that the program tests to split data with itself
    assert( newSplitFeatureIdx != targetIdx );

    vector<size_t> newSampleIcs_left(0);
    vector<size_t> newSampleIcs_right = sampleIcs;
    num_t newSplitValue;
    set<num_t> newSplitValues_left;
    set<num_t> newSplitValues_right;
    num_t newSplitFitness = 0.0;

    if ( isFeatureNumerical ) {

      newSplitFitness = treeData->numericalFeatureSplit(targetIdx,
							newSplitFeatureIdx,
							parameters->nodeSize,
							newSampleIcs_left,
							newSampleIcs_right,
							newSplitValue);
    } else {

      newSplitFitness = treeData->categoricalFeatureSplit(targetIdx,
							  newSplitFeatureIdx,
							  parameters->nodeSize,
							  newSampleIcs_left,
							  newSampleIcs_right,
							  newSplitValues_left,
							  newSplitValues_right);
    }

    if( newSplitFitness > splitFitness &&
	newSampleIcs_left.size() >= parameters->nodeSize &&
	newSampleIcs_right.size() >= parameters->nodeSize ) {
      
      splitFitness = newSplitFitness;
      splitFeatureIdx = newSplitFeatureIdx;
      splitValue = newSplitValue;
      splitValues_left = newSplitValues_left;
      splitValues_right = newSplitValues_right;
      sampleIcs_left = newSampleIcs_left;
      sampleIcs_right = newSampleIcs_right;
    }    

  }
  
  // If none of the splitter candidates worked as a splitter
  if ( fabs(splitFitness) < datadefs::EPS ) {
    return(false);
  }

  // Get the fraction of train samples that were sent to the left child
  num_t leftFraction = 1.0 * sampleIcs_left.size() / ( sampleIcs_left.size() + sampleIcs_right.size() );
  
  if ( treeData->isFeatureNumerical(splitFeatureIdx) ) {

    this->setSplitter(splitFeatureIdx,treeData->getFeatureName(splitFeatureIdx),leftFraction,splitValue);

  } else {
    
    set<string> rawSplitValues_left,rawSplitValues_right;

    for ( set<num_t>::const_iterator it(splitValues_left.begin()); it != splitValues_left.end(); ++it ) {
      rawSplitValues_left.insert( treeData->getRawFeatureData(splitFeatureIdx,*it) );
    }

    for ( set<num_t>::const_iterator it(splitValues_right.begin()); it != splitValues_right.end(); ++it ) {
      rawSplitValues_right.insert( treeData->getRawFeatureData(splitFeatureIdx,*it) );
    }

    this->setSplitter(splitFeatureIdx,treeData->getFeatureName(splitFeatureIdx),leftFraction,rawSplitValues_left,rawSplitValues_right);

  }

  return(true);

}


  
