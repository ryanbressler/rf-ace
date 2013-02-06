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
		       const unordered_set<string>& leftSplitValues, 
		       Node& leftChild,
		       Node& rightChild) {

  if ( this->hasChildren() ) {
    cerr << "Node::setSplitter() -- cannot set a splitter to a node twice!" << endl;
    exit(1);
  }

  splitter_.name = splitterName;
  splitter_.type = Feature::Type::CAT;
  splitter_.leftValues = leftSplitValues;

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
      } else {
	return( this->rightChild()->percolate(testData,sampleIdx,scrambleFeatureIdx) );
      }
    }
    
  } else {
    
    if ( testData->feature(featureIdx)->hasHash(sampleIdx,splitter_.hashValue) ) {
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
	     << ",LVALUES=" << splitter_.leftLeqValue; 
      //<< ",RVALUES=" << splitter_.leftLeqValue;
    } else if ( splitter_.type == Feature::Type::CAT ) {
      toFile << ",SPLITTERTYPE=CATEGORICAL" 
	     << ",LVALUES=" << "\""; utils::write(toFile,splitter_.leftValues.begin(),splitter_.leftValues.end(),':'); toFile << "\"";
      //<< ",RVALUES=" << "\""; utils::write(toFile,splitter_.rightValues.begin(),splitter_.rightValues.end(),':'); toFile << "\"";
    } else {
      toFile << ",SPLITTERTYPE=TEXTUAL"
	     << ",LVALUES=" << splitter_.hashValue;
      //<< ",RVALUES=" << splitter_.hashValue;
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
			      vector<Node>& children,
			      SplitCache& splitCache) {

  if ( false ) {
    cout << "REC " << this << endl;
    cout << " " << treeData << endl;
    cout << " " << forestOptions << endl;
    cout << " " << nLeaves << endl;
  }

  splitCache.trainData = treeData->getFeatureData(targetIdx,sampleIcs);
  splitCache.nSamples = sampleIcs.size();

  if ( predictionFunctionType == MEAN ) {
    num_t trainPrediction = math::mean(splitCache.trainData);
    string rawTrainPrediction = utils::num2str(trainPrediction);
    this->setTrainPrediction( trainPrediction, rawTrainPrediction );
  } else if ( predictionFunctionType == MODE ) {
    num_t trainPrediction = math::mode<num_t>(splitCache.trainData);
    string rawTrainPrediction = treeData->getRawFeatureData(targetIdx,trainPrediction);
    this->setTrainPrediction( trainPrediction, rawTrainPrediction );
  } else if ( predictionFunctionType == GAMMA ) {
    num_t trainPrediction = math::gamma(splitCache.trainData, treeData->feature(targetIdx)->nCategories() );
    string rawTrainPrediction = utils::num2str(trainPrediction);
    this->setTrainPrediction( trainPrediction, rawTrainPrediction );
  } else {
    cerr << "Node::recursiveNodeSplit() -- unknown prediction function!" << endl;
    exit(1);
  }

  assert( !datadefs::isNAN(trainPrediction_) );
  assert( !datadefs::isNAN_STR(rawTrainPrediction_) );

  assert( *nLeaves <= forestOptions->nMaxLeaves );

  if ( splitCache.nSamples < 2 * forestOptions->nodeSize || *nLeaves == forestOptions->nMaxLeaves || childIdx + 1 >= children.size() ) {
    return;
  }

  splitCache.featureSampleIcs.clear();

  if ( forestOptions->isRandomSplit ) {

    splitCache.featureSampleIcs.resize(forestOptions->mTry);

    for ( size_t i = 0; i < forestOptions->mTry; ++i ) {
      splitCache.featureSampleIcs[i] = pmf->sample(random); //icdf( random->uniform() );
    }

    if ( forestOptions->useContrasts ) {
      for ( size_t i = 0; i < forestOptions->mTry; ++i ) {
	
	// If the sampled feature is a contrast... 
	if ( ! treeData->feature(splitCache.featureSampleIcs[i])->isTextual() && random->uniform() < forestOptions->contrastFraction ) { // p% sampling rate
	  
	  // Contrast features in Treedata are indexed with an offset of the number of features: nFeatures
	  splitCache.featureSampleIcs[i] += treeData->nFeatures();
	}
      }
    } 
  } else {

    splitCache.featureSampleIcs = utils::range(treeData->nFeatures());

    splitCache.featureSampleIcs.erase(splitCache.featureSampleIcs.begin()+targetIdx);
  }
  
  splitCache.sampleIcs_left.clear();
  splitCache.sampleIcs_right.clear();
  splitCache.sampleIcs_missing.clear();

  bool foundSplit = this->regularSplitterSeek(treeData,
					      targetIdx,
					      forestOptions,
					      random,
					      sampleIcs,
					      childIdx,
					      children,
					      splitCache);
        
  if ( !foundSplit ) {
    return;
  }
  
  if ( minDistToRoot[splitCache.splitFeatureIdx] > treeDepth ) {
    minDistToRoot[splitCache.splitFeatureIdx] = treeDepth;
  }
  
  featuresInTree.insert(splitCache.splitFeatureIdx);
  
  *nLeaves += 1;

  vector<size_t> sampleIcs_left = splitCache.sampleIcs_left;
  vector<size_t> sampleIcs_right = splitCache.sampleIcs_right;
  vector<size_t> sampleIcs_missing = splitCache.sampleIcs_missing;
  
  // Left child recursive split
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
					children,
					splitCache);


  // Right child recursive split
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
					 children,
					 splitCache);
  

  
  // OPTIONAL: Missing child recursive split
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
					     children,
					     splitCache); 
  }
  
}

bool Node::regularSplitterSeek(Treedata* treeData,
			       const size_t targetIdx,
			       const ForestOptions* forestOptions,
			       distributions::Random* random,
			       const vector<size_t>& sampleIcs,
			       size_t& childIdx,
			       vector<Node>& children,
			       SplitCache& splitCache) {
  
  // This many features will be tested for splitting the data
  size_t nFeaturesForSplit = splitCache.featureSampleIcs.size();
  
  // Initialize split fitness to lowest possible value
  splitCache.splitFitness = 0.0;

  // Loop through candidate splitters
  for ( size_t i = 0; i < nFeaturesForSplit; ++i ) {
    
    // Get split feature index
    splitCache.newSplitFeatureIdx = splitCache.featureSampleIcs[i];

    // We don't want that the program tests to split data with itself
    assert( splitCache.newSplitFeatureIdx != targetIdx );

    // Reset the splitCache
    splitCache.newSampleIcs_left.clear();
    splitCache.newSampleIcs_right = sampleIcs;
    splitCache.newSampleIcs_missing.clear();
    treeData->separateMissingSamples(splitCache.newSplitFeatureIdx,splitCache.newSampleIcs_right,splitCache.newSampleIcs_missing);
    splitCache.newSplitValue = datadefs::NUM_NAN;
    splitCache.newSplitValues_left.clear();
    splitCache.newHashIdx = 0;
    splitCache.newSplitFitness = 0.0;

    const Feature* newSplitFeature = treeData->feature(splitCache.newSplitFeatureIdx);

    if ( newSplitFeature->isNumerical() ) {

      splitCache.newSplitFitness = treeData->numericalFeatureSplit(targetIdx,
								   splitCache.newSplitFeatureIdx,
								   forestOptions->nodeSize,
								   splitCache.newSampleIcs_left,
								   splitCache.newSampleIcs_right,
								   splitCache.newSplitValue);

    } else if ( newSplitFeature->isCategorical() ) {
      
      unordered_set<num_t> uniqueCats(sampleIcs.size());
      
      for ( size_t i = 0; i < splitCache.newSampleIcs_right.size(); ++i ) {
	uniqueCats.insert(treeData->feature(splitCache.newSplitFeatureIdx)->data[splitCache.newSampleIcs_right[i]]);
      }
      
      vector<num_t> catOrder(uniqueCats.size());
      size_t iter = 0;
      for ( unordered_set<num_t>::const_iterator it(uniqueCats.begin()); it != uniqueCats.end(); ++it ) {
	catOrder[iter] = *it;
	++iter;
      }
      
      utils::permute(catOrder,random);
      
      splitCache.newSplitFitness = treeData->categoricalFeatureSplit(targetIdx,
								     splitCache.newSplitFeatureIdx,
								     catOrder,
								     forestOptions->nodeSize,
								     splitCache.newSampleIcs_left,
								     splitCache.newSampleIcs_right,
								     splitCache.newSplitValues_left);

    } else if ( newSplitFeature->isTextual() && splitCache.newSampleIcs_right.size() > 0 ) {

      // Choose random sample
      size_t sampleIdx = splitCache.newSampleIcs_right[ random->integer() % splitCache.newSampleIcs_right.size() ];

      // Choose random hash from the randomly selected sample
      splitCache.newHashIdx = newSplitFeature->getHash(sampleIdx,random->integer());

      splitCache.newSplitFitness = treeData->textualFeatureSplit(targetIdx,
								 splitCache.newSplitFeatureIdx,
								 splitCache.newHashIdx,
								 forestOptions->nodeSize,
								 splitCache.newSampleIcs_left,
								 splitCache.newSampleIcs_right);

    }

    if( splitCache.newSplitFitness > splitCache.splitFitness &&
	splitCache.newSampleIcs_left.size() >= forestOptions->nodeSize &&
	splitCache.newSampleIcs_right.size() >= forestOptions->nodeSize ) {
      
      splitCache.splitFitness      = splitCache.newSplitFitness;
      splitCache.splitFeatureIdx   = splitCache.newSplitFeatureIdx;
      splitCache.splitValue        = splitCache.newSplitValue;
      splitCache.splitValues_left  = splitCache.newSplitValues_left;
      splitCache.hashIdx           = splitCache.newHashIdx;
      splitCache.sampleIcs_left    = splitCache.newSampleIcs_left;
      splitCache.sampleIcs_right   = splitCache.newSampleIcs_right;
      splitCache.sampleIcs_missing = splitCache.newSampleIcs_missing;
    }    

  }
  
  // If none of the splitter candidates worked as a splitter
  if ( fabs(splitCache.splitFitness) < datadefs::EPS ) {
    return(false);
  } 

  const Feature* splitFeature = treeData->feature(splitCache.splitFeatureIdx);

  if ( splitFeature->isNumerical() ) {

    this->setSplitter(splitFeature->name(),splitCache.splitValue,children[childIdx],children[childIdx+1]);

  } else if ( splitFeature->isCategorical() ) {
    
    unordered_set<string> rawSplitValues_left;

    for ( unordered_set<num_t>::const_iterator it(splitCache.splitValues_left.begin()); it != splitCache.splitValues_left.end(); ++it ) {
      rawSplitValues_left.insert( treeData->getRawFeatureData(splitCache.splitFeatureIdx,*it) );
    }


    this->setSplitter(splitFeature->name(),rawSplitValues_left,children[childIdx],children[childIdx+1]);

  } else if ( splitFeature->isTextual() ) {
    
    this->setSplitter(splitFeature->name(),splitCache.hashIdx,children[childIdx],children[childIdx+1]);

  }

  childIdx += 2;

  if ( ! forestOptions->noNABranching && splitCache.sampleIcs_missing.size() > 0 && childIdx < children.size() ) { 
    missingChild_ = &children[childIdx++];
  }

  return(true);

}


  
