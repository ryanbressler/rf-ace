#include<iostream>
#include<cassert>
#include<iomanip>

#include "node.hpp"
#include "utils.hpp"
#include "math.hpp"

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

void Node::GrowInstructions::validate() const {

  assert( maxNodesToStop > 0 );
  assert( minNodeSizeToStop > 0 );
  assert( nFeaturesForSplit > 0 );
  assert( featureIcs.size() > 0 );
  
  if ( isRandomSplit ) {
    assert( nFeaturesForSplit <= featureIcs.size() );
  } else {
    assert( nFeaturesForSplit == featureIcs.size() );
  }

  if ( sampleWithReplacement ) {
    assert( 0.0 < sampleSizeFraction );
  } else {
    assert( 0.0 < sampleSizeFraction && sampleSizeFraction <= 1.0 );
  }

}


/**
 * Deletes child nodes, which will cascade all the way to the leaf nodes 
 */
void Node::deleteTree() {
  
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

  toFile << "NODE=" << traversal << ",PRED=" << trainPrediction_;
  
  if ( this->hasChildren() ) {
    
    string splitterType = splitter_.isNumerical ? "NUMERICAL" : "CATEGORICAL";

    toFile << ",SPLITTER=" << "\"" << splitter_.name << "\""
	   << ",SPLITTERTYPE=" << splitterType
	   << ",LFRACTION=" << splitter_.leftFraction; 

    if ( splitter_.isNumerical ) {
      toFile << ",LVALUES=" << splitter_.leftLeqValue
	     << ",RVALUES=" << splitter_.leftLeqValue << endl;
    } else {

      string leftValues = utils::join(splitter_.leftValues.begin(),splitter_.leftValues.end(),':');
      string rightValues = utils::join(splitter_.rightValues.begin(),splitter_.rightValues.end(),':');
      
      toFile << ",LVALUES=" << "\"" << leftValues << "\""
	     << ",RVALUES=" << "\"" << rightValues << "\"" << endl;
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
                              const vector<size_t>& sampleIcs,
                              const GrowInstructions& GI,
                              set<size_t>& featuresInTree,
                              size_t* nNodes) {

  size_t nSamples = sampleIcs.size();

  assert( *nNodes <= GI.maxNodesToStop );

  if ( nSamples < 2 * GI.minNodeSizeToStop || *nNodes >= GI.maxNodesToStop ) {
    
    vector<num_t> trainData = treeData->getFeatureData(targetIdx,sampleIcs);

    if ( GI.predictionFunctionType == MEAN ) {
      this->setTrainPrediction( math::mean(trainData) );
    } else if ( GI.predictionFunctionType == MODE ) {
      num_t trainPrediction = math::mode<num_t>(trainData);
      string rawTrainPrediction = treeData->getRawFeatureData(targetIdx,trainPrediction);
      this->setTrainPrediction( trainPrediction, rawTrainPrediction );
    } else if ( GI.predictionFunctionType == GAMMA ) {
      this->setTrainPrediction( math::gamma(trainData,GI.numClasses) );
    } else {
      cerr << "Node::recursiveNodeSplit() -- unknown prediction function!" << endl;
      exit(1);
    }

    assert( !datadefs::isNAN(trainPrediction_) );

    return;
  }

  vector<size_t> featureSampleIcs = GI.featureIcs;

  if(GI.isRandomSplit) {

    // In-place permute of feature indices
    treeData->permute<size_t>(featureSampleIcs);

    // Take only the first ones
    featureSampleIcs.resize(GI.nFeaturesForSplit);

  }

  // With 1% sampling rate assign contrasts
  if(GI.useContrasts) {
    for(size_t i = 0; i < GI.nFeaturesForSplit; ++i) {
      
      // If the sampled feature is a contrast... 
      if( treeData->getRandomUnif() < 0.01 ) { // %1 sampling rate
	
	featureSampleIcs[i] += treeData->nFeatures();
      }
    }
  } 

  //datadefs::print<size_t>(featureSampleIcs);
  
  assert( featureSampleIcs.size() == GI.nFeaturesForSplit );
  
  vector<size_t> sampleIcs_left,sampleIcs_right;

  size_t splitFeatureIdx;
  num_t splitFitness;
  
  bool foundSplit = this->regularSplitterSeek(treeData,
					      targetIdx,
					      sampleIcs,
					      featureSampleIcs,
					      GI,
					      splitFeatureIdx,
					      sampleIcs_left,
					      sampleIcs_right,
					      splitFitness);
        
  if ( !foundSplit ) {

    vector<num_t> trainData = treeData->getFeatureData(targetIdx,sampleIcs);
    
    if ( GI.predictionFunctionType == MEAN ) {
      this->setTrainPrediction( math::mean(trainData) );
    } else if ( GI.predictionFunctionType == MODE ) {
      num_t trainPrediction = math::mode<num_t>(trainData);
      string rawTrainPrediction = treeData->getRawFeatureData(targetIdx,trainPrediction);
      this->setTrainPrediction( trainPrediction, rawTrainPrediction );
    } else if ( GI.predictionFunctionType == GAMMA ) {
      this->setTrainPrediction( math::gamma(trainData,GI.numClasses) );
    } else {
      cerr << "Node::recursiveNodeSplit() -- unknown prediction function!" << endl;
      exit(1);
    }
    
    assert( !datadefs::isNAN(trainPrediction_) );

    return;
  }
  
  featuresInTree.insert(splitFeatureIdx);
  *nNodes += 2;
  
  leftChild_->recursiveNodeSplit(treeData,targetIdx,sampleIcs_left,GI,featuresInTree,nNodes);
  rightChild_->recursiveNodeSplit(treeData,targetIdx,sampleIcs_right,GI,featuresInTree,nNodes);
  
}

bool Node::regularSplitterSeek(Treedata* treeData,
			       const size_t targetIdx,
			       const vector<size_t>& sampleIcs,
			       const vector<size_t>& featureSampleIcs,
			       const GrowInstructions& GI,
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
							GI.minNodeSizeToStop,
							newSampleIcs_left,
							newSampleIcs_right,
							newSplitValue);
    } else {

      newSplitFitness = treeData->categoricalFeatureSplit(targetIdx,
							  newSplitFeatureIdx,
							  GI.minNodeSizeToStop,
							  newSampleIcs_left,
							  newSampleIcs_right,
							  newSplitValues_left,
							  newSplitValues_right);
    }

    if( newSplitFitness > splitFitness &&
	newSampleIcs_left.size() >= GI.minNodeSizeToStop &&
	newSampleIcs_right.size() >= GI.minNodeSizeToStop ) {
      
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


  
