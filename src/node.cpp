#include<iostream>
#include<cassert>
#include<iomanip>

#include "node.hpp"
#include "datadefs.hpp"
//#include "utils.hpp"
#include "math.hpp"

using namespace std;
using datadefs::num_t;

Node::Node():
  trainPrediction_(datadefs::NUM_NAN),
  rawTrainPrediction_(datadefs::STR_NAN),
  leftChild_(NULL),
  rightChild_(NULL),
  missingChild_(NULL) {
}

Node::~Node() { }

// !! Documentation: consider combining with the documentation in the header
// !! file, fleshing it out a bit. Ideally, implementation notes should fall
// !! here; notes on the abstraction should fall in the header file.
void Node::setSplitter(const string& splitterName,  
		       const num_t splitLeftLeqValue, 
		       Node& leftChild, 
		       Node& rightChild) {
  
  if ( this->hasChildren() ) {
    cerr << "Cannot set a splitter to a node twice!" << endl;
    exit(1);
  }

  splitter_.name = splitterName;
  splitter_.type = Feature::Type::NUM;
  splitter_.leftLeqValue = splitLeftLeqValue;

  leftChild_ = &leftChild;
  rightChild_ = &rightChild;

}

// !! Documentation: consider combining with the documentation in the header
// !! file, fleshing it out a bit. Ideally, implementation notes should fall
// !! here; notes on the abstraction should fall in the header file.
void Node::setSplitter(const string& splitterName,  
		       const set<string>& leftSplitValues, 
		       const set<string>& rightSplitValues,
		       Node& leftChild,
		       Node& rightChild) {

  if ( this->hasChildren() ) {
    cerr << "Node::setSplitter() -- cannot set a splitter to a node twice!" << endl;
    exit(1);
  }

  splitter_.name = splitterName;
  splitter_.type = Feature::Type::CAT;
  splitter_.leftValues = leftSplitValues;
  splitter_.rightValues = rightSplitValues;

  leftChild_ = &leftChild;
  rightChild_ = &rightChild;

}

void Node::setSplitter(const string& splitterName,
		       const uint32_t hashIdx,
		       Node& leftChild,
		       Node& rightChild) {

  if ( this->hasChildren() ) {
    cerr << "Node::setSplitter() -- cannot set a splitter to a node twice!" << endl;
    exit(1);
  }

  splitter_.name = splitterName;
  splitter_.type = Feature::Type::TXT;
  splitter_.hashValue = hashIdx;

  leftChild_ = &leftChild;
  rightChild_ = &rightChild;

}

void Node::setMissingChild(Node& missingChild) {
  missingChild_ = &missingChild;
}

const Node* Node::percolate(Treedata* testData, const size_t sampleIdx, const size_t scrambleFeatureIdx) const {
  
  if ( !this->hasChildren() ) { return( this ); }
  
  size_t featureIdx = testData->getFeatureIdx(splitter_.name);
  
  if ( featureIdx == testData->end() ) { return( this ); }
  
  if ( splitter_.type == Feature::Type::NUM ) {
    num_t data;
    if ( scrambleFeatureIdx != featureIdx ) {
      data = testData->getFeatureData(featureIdx,sampleIdx);
    } else {
      cerr << "Randomized prediction is not available!" << endl;
      exit(1);
      // data = testData->getFeatureData(featureIdx, random_->integer() % testData->nSamples() );
    }
    if ( datadefs::isNAN(data) ) { 
      if ( this->missingChild() ) {
	return( this->missingChild()->percolate(testData,sampleIdx,scrambleFeatureIdx) );
      } else {
	return( this ); 
      }
    } else {
      return( data <= splitter_.leftLeqValue ? 
	      this->leftChild()->percolate(testData,sampleIdx,scrambleFeatureIdx) : 
	      this->rightChild()->percolate(testData,sampleIdx,scrambleFeatureIdx) );
    }
  } else if ( splitter_.type == Feature::Type::CAT ){
    string data;
    if ( scrambleFeatureIdx != featureIdx ) {
      data = testData->getRawFeatureData(featureIdx,sampleIdx);
    } else {
      cerr << "randomized prediction is not available!" << endl;
      exit(1);
      // data = testData->getRawFeatureData(featureIdx, random_->integer() % testData->nSamples() );
    }
    if ( datadefs::isNAN_STR(data) ) { 
      if ( this->missingChild() ) {
	return( this->missingChild()->percolate(testData,sampleIdx,scrambleFeatureIdx) );
      } else {
	return( this );
      }
    } else { 
      
      // Return left child if splits left
      if ( splitter_.leftValues.find(data) != splitter_.leftValues.end() ) {
	return( this->leftChild()->percolate(testData,sampleIdx,scrambleFeatureIdx) );
      }
      // Return right child if splits right
      if ( splitter_.rightValues.find(data) != splitter_.rightValues.end() ) {
	return( this->rightChild()->percolate(testData,sampleIdx,scrambleFeatureIdx) );
      }
      
      // Else return this
      return( this );
    }
    
  } else {
    
    if ( testData->hasHash(featureIdx,sampleIdx,splitter_.hashValue) ) {
      return( this->leftChild()->percolate(testData,sampleIdx,scrambleFeatureIdx) );
    } else {
      return( this->rightChild()->percolate(testData,sampleIdx,scrambleFeatureIdx) );
    }
  }
  
}

Node* Node::leftChild() const {
  return( leftChild_ );
}


Node* Node::rightChild() const {
  return( rightChild_ );
}

Node* Node::missingChild() const {
  return( missingChild_ );
}

/**
 * Recursively prints a tree to a stream (file)
 */
void Node::print(string& traversal, ofstream& toFile) {

  toFile << "NODE=" << traversal << ",PRED=" << rawTrainPrediction_;
  
  if ( this->hasChildren() ) {
    
    toFile << ",SPLITTER=" << "\"" << splitter_.name << "\"";

    if (splitter_.type == Feature::Type::NUM ) {
      toFile << ",SPLITTERTYPE=NUMERICAL"
	     << ",LVALUES=" << splitter_.leftLeqValue 
	     << ",RVALUES=" << splitter_.leftLeqValue;
    } else if ( splitter_.type == Feature::Type::CAT ) {
      toFile << ",SPLITTERTYPE=CATEGORICAL" 
	     << ",LVALUES=" << "\""; utils::write(toFile,splitter_.leftValues.begin(),splitter_.leftValues.end(),':'); toFile << "\""
	     << ",RVALUES=" << "\""; utils::write(toFile,splitter_.rightValues.begin(),splitter_.rightValues.end(),':'); toFile << "\"";
    } else {
      toFile << ",SPLITTERTYPE=TEXTUAL"
	     << ",LVALUES=" << splitter_.hashValue
	     << ",RVALUES=" << splitter_.hashValue;
    }

    if ( this->missingChild() ) { toFile << ",M=M" << endl; } else { toFile << endl; }

    string traversalLeft = traversal;
    traversalLeft.append("L");
    string traversalRight = traversal;
    traversalRight.append("R");
    
    this->leftChild()->print(traversalLeft,toFile);
    this->rightChild()->print(traversalRight,toFile);

    // Optional third branch
    if ( this->missingChild() ) {
      string traversalMissing = traversal;
      traversalMissing.append("M");
      this->missingChild()->print(traversalMissing,toFile);
    }
    
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
num_t Node::getTrainPrediction() const {
  assert( !datadefs::isNAN(trainPrediction_) );
  return( trainPrediction_ );
}

string Node::getRawTrainPrediction() const {
  assert( !datadefs::isNAN_STR(rawTrainPrediction_) );
  return( rawTrainPrediction_ );
}

void Node::recursiveNodeSplit(Treedata* treeData,
			      const size_t targetIdx,
			      const ForestOptions* forestOptions,
			      distributions::Random* random,
			      const PredictionFunctionType& predictionFunctionType,
			      const distributions::PMF* pmf,
			      const vector<size_t>& sampleIcs,
			      const size_t treeDepth,
			      set<size_t>& featuresInTree,
			      vector<size_t>& minDistToRoot,
			      size_t* nLeaves,
			      size_t& childIdx,
			      vector<Node>& children) {

  if ( false ) {
    cout << "REC " << this << endl;
    cout << " " << treeData << endl;
    cout << " " << forestOptions << endl;
    cout << " " << nLeaves << endl;
  }

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

  assert( !datadefs::isNAN(trainPrediction_) );
  assert( !datadefs::isNAN_STR(rawTrainPrediction_) );

  size_t nSamples = sampleIcs.size();

  assert( *nLeaves <= forestOptions->nMaxLeaves );

  if ( nSamples < 2 * forestOptions->nodeSize || *nLeaves == forestOptions->nMaxLeaves || childIdx + 3 >= children.size() ) {
    return;
  }

  vector<size_t> featureSampleIcs;

  if ( forestOptions->isRandomSplit ) {

    featureSampleIcs.resize(forestOptions->mTry);

    for ( size_t i = 0; i < forestOptions->mTry; ++i ) {
      featureSampleIcs[i] = pmf->icdf( random->uniform() );
    }

    if ( forestOptions->useContrasts ) {
      for ( size_t i = 0; i < forestOptions->mTry; ++i ) {
	
	// If the sampled feature is a contrast... 
	if ( ! treeData->isFeatureTextual(featureSampleIcs[i]) && random->uniform() < forestOptions->contrastFraction ) { // p% sampling rate
	  
	  // Contrast features in Treedata are indexed with an offset of the number of features: nFeatures
	  featureSampleIcs[i] += treeData->nFeatures();
	}
      }
    } 
  } else {

    featureSampleIcs = utils::range(treeData->nFeatures());

    featureSampleIcs.erase(featureSampleIcs.begin()+targetIdx);
  }
  // assert( featureSampleIcs.size() == forestOptions->mTry );
  
  vector<size_t> sampleIcs_left,sampleIcs_right,sampleIcs_missing;

  size_t splitFeatureIdx;
  num_t splitFitness;

  bool foundSplit = this->regularSplitterSeek(treeData,
					      targetIdx,
					      forestOptions,
					      random,
					      sampleIcs,
					      featureSampleIcs,
					      splitFeatureIdx,
					      sampleIcs_left,
					      sampleIcs_right,
					      sampleIcs_missing,
					      splitFitness,
					      childIdx,
					      children);
        
  if ( !foundSplit ) {
    return;
  }

  if ( minDistToRoot[splitFeatureIdx] > treeDepth ) {
    minDistToRoot[splitFeatureIdx] = treeDepth;
  }

  featuresInTree.insert(splitFeatureIdx);

  *nLeaves += 1;
  
  this->leftChild()->recursiveNodeSplit(treeData,
					targetIdx,
					forestOptions,
					random,
					predictionFunctionType,
					pmf,
					sampleIcs_left,
					treeDepth+1,
					featuresInTree,
					minDistToRoot,
					nLeaves,
					childIdx,
					children);

  this->rightChild()->recursiveNodeSplit(treeData,
					 targetIdx,
					 forestOptions,
					 random,
					 predictionFunctionType,
					 pmf,
					 sampleIcs_right,
					 treeDepth+1,
					 featuresInTree,
					 minDistToRoot,
					 nLeaves,
					 childIdx,
					 children);

  if ( this->missingChild() ) {

    assert( sampleIcs_missing.size() > 0 );

    *nLeaves += 1;
    
    this->missingChild()->recursiveNodeSplit(treeData,
					     targetIdx,
					     forestOptions,
					     random,
					     predictionFunctionType,
					     pmf,
					     sampleIcs_missing,
					     treeDepth+1,
					     featuresInTree,
					     minDistToRoot,
					     nLeaves,
					     childIdx,
					     children);
    
  }
  
}

bool Node::regularSplitterSeek(Treedata* treeData,
			       const size_t targetIdx,
			       const ForestOptions* forestOptions,
			       distributions::Random* random,
			       const vector<size_t>& sampleIcs,
			       const vector<size_t>& featureSampleIcs,
			       size_t& splitFeatureIdx,
			       vector<size_t>& sampleIcs_left,
			       vector<size_t>& sampleIcs_right,
			       vector<size_t>& sampleIcs_missing,
			       num_t& splitFitness,
			       size_t& childIdx,
			       vector<Node>& children) {
  
  // This many features will be tested for splitting the data
  size_t nFeaturesForSplit = featureSampleIcs.size();
  
  num_t splitValue = datadefs::NUM_NAN;
  set<num_t> splitValues_left; 
  set<num_t> splitValues_right;
  uint32_t hashIdx = 0;

  // Initialize split fitness to lowest possible value
  splitFitness = 0.0;

  // Loop through candidate splitters
  for ( size_t i = 0; i < nFeaturesForSplit; ++i ) {
    
    // Get split feature index
    size_t newSplitFeatureIdx = featureSampleIcs[i];

    // We don't want that the program tests to split data with itself
    assert( newSplitFeatureIdx != targetIdx );

    vector<size_t> newSampleIcs_left(0);
    vector<size_t> newSampleIcs_right = sampleIcs;
    vector<size_t> newSampleIcs_missing(0);
    num_t newSplitValue;
    set<num_t> newSplitValues_left;
    set<num_t> newSplitValues_right;
    uint32_t newHashIdx = 0;
    num_t newSplitFitness = 0.0;

    if ( treeData->isFeatureNumerical(newSplitFeatureIdx) ) {

      newSplitFitness = treeData->numericalFeatureSplit(targetIdx,
							newSplitFeatureIdx,
							forestOptions->nodeSize,
							newSampleIcs_left,
							newSampleIcs_right,
							newSampleIcs_missing,
							newSplitValue);

    } else if ( treeData->isFeatureCategorical(newSplitFeatureIdx) ) {

      newSplitFitness = treeData->categoricalFeatureSplit(targetIdx,
							  newSplitFeatureIdx,
							  forestOptions->nodeSize,
							  newSampleIcs_left,
							  newSampleIcs_right,
							  newSampleIcs_missing,
							  newSplitValues_left,
							  newSplitValues_right);
    } else if ( treeData->isFeatureTextual(newSplitFeatureIdx) && newSampleIcs_right.size() > 0 ) {

      // Choose random sample
      size_t sampleIdx = random->integer() % newSampleIcs_right.size();

      // Choose random hash from the randomly selected sample
      newHashIdx = treeData->getHash(newSplitFeatureIdx,sampleIdx,random->integer());

      newSplitFitness = treeData->textualFeatureSplit(targetIdx,
						      newSplitFeatureIdx,
						      newHashIdx,
						      forestOptions->nodeSize,
						      newSampleIcs_left,
						      newSampleIcs_right,
						      newSampleIcs_missing);

    }

    if( newSplitFitness > splitFitness &&
	newSampleIcs_left.size() >= forestOptions->nodeSize &&
	newSampleIcs_right.size() >= forestOptions->nodeSize ) {
      
      splitFitness = newSplitFitness;
      splitFeatureIdx = newSplitFeatureIdx;
      splitValue = newSplitValue;
      splitValues_left = newSplitValues_left;
      splitValues_right = newSplitValues_right;
      hashIdx = newHashIdx;
      sampleIcs_left = newSampleIcs_left;
      sampleIcs_right = newSampleIcs_right;
      sampleIcs_missing = newSampleIcs_missing;
    }    

  }
  
  // If none of the splitter candidates worked as a splitter
  if ( fabs(splitFitness) < datadefs::EPS ) {
    return(false);
  } 

  if ( treeData->isFeatureNumerical(splitFeatureIdx) ) {

    this->setSplitter(treeData->getFeatureName(splitFeatureIdx),splitValue,children[childIdx],children[childIdx+1]);

  } else if ( treeData->isFeatureCategorical(splitFeatureIdx) ) {
    
    set<string> rawSplitValues_left, rawSplitValues_right;

    for ( set<num_t>::const_iterator it(splitValues_left.begin()); it != splitValues_left.end(); ++it ) {
      rawSplitValues_left.insert( treeData->getRawFeatureData(splitFeatureIdx,*it) );
    }

    for ( set<num_t>::const_iterator it(splitValues_right.begin()); it != splitValues_right.end(); ++it ) {
      rawSplitValues_right.insert( treeData->getRawFeatureData(splitFeatureIdx,*it) );
    }

    this->setSplitter(treeData->getFeatureName(splitFeatureIdx),rawSplitValues_left,rawSplitValues_right,children[childIdx],children[childIdx+1]);

  } else if ( treeData->isFeatureTextual(splitFeatureIdx) ) {
    
    this->setSplitter(treeData->getFeatureName(splitFeatureIdx),hashIdx,children[childIdx],children[childIdx+1]);

  }

  childIdx += 2;

  if ( ! forestOptions->noNABranching && sampleIcs_missing.size() > 0 ) { 
    missingChild_ = &children[childIdx++];
  }

  return(true);

}


  
