#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "stochasticforest.hpp"
#include "datadefs.hpp"
#include "argparse.hpp"
#include "utils.hpp"

StochasticForest::StochasticForest(Treedata* treeData, const string& targetName, const Parameters& parameters):
  treeData_(treeData),
  parameters_(parameters),
  targetName_(targetName),
  rootNodes_(parameters_.nTrees),
  importanceValues_(treeData_->nFeatures(),0.0),
  contrastImportanceValues_(treeData_->nFeatures(),0.0),
  oobError_(datadefs::NUM_NAN),
  oobMatrix_(parameters_.nTrees) {
  
  featuresInForest_.clear();

  size_t targetIdx = treeData_->getFeatureIdx(targetName_);
  targetSupport_ = treeData_->categories(targetIdx);

  if ( parameters_.model == RF ) {
    this->learnRF();
  } else {
    this->learnGBT();
  }

}

StochasticForest::StochasticForest(Treedata* treeData, const string& forestFile):
  treeData_(treeData) {

  ifstream forestStream( forestFile.c_str() );
  assert(forestStream.good());
  
  string newLine("");
  getline(forestStream,newLine);

  map<string,string> forestSetup = utils::parse(newLine,',','=','"');

  //assert( forestSetup["FOREST"] == "GBT" );
  if ( forestSetup["FOREST"] == "GBT" ) {
    parameters_.model = GBT;
  } else if ( forestSetup["FOREST"] == "RF" ) {
    parameters_.model = RF;
  } else {
    cerr << "Unknown forest type: " << forestSetup["FOREST"] << endl;
    exit(1);
  }

  assert( forestSetup.find("NTREES") != forestSetup.end() );
  parameters_.nTrees = static_cast<size_t>( utils::str2int(forestSetup["NTREES"]) );
  rootNodes_.resize(parameters_.nTrees);

  assert( forestSetup.find("TARGET") != forestSetup.end() );
  targetName_ = forestSetup["TARGET"];

  assert( forestSetup.find("CATEGORIES") != forestSetup.end() );
  targetSupport_ = utils::split(forestSetup["CATEGORIES"],',');

  assert( forestSetup.find("SHRINKAGE") != forestSetup.end() );
  parameters_.shrinkage = datadefs::str2num(forestSetup["SHRINKAGE"]);
  
  assert( forestStream.good() );
 
  int treeIdx = -1;

  vector<map<string,Node*> > forestMap( parameters_.nTrees );

  // Read the forest
  while ( getline(forestStream,newLine) ) {

    // remove trailing end-of-line characters
    newLine = utils::chomp(newLine);

    // Empty lines will be discarded
    if ( newLine == "" ) {
      continue;
    }

    if ( newLine.compare(0,5,"TREE=") == 0 ) { 
      ++treeIdx;
      //cout << treeIdx << " vs " << newLine << " => " << utils::str2int(utils::split(newLine,'=')[1]) << endl;
      assert( treeIdx == utils::str2int(utils::split(newLine,'=')[1]) );
      continue;
    }
        
    //cout << newLine << endl;

    map<string,string> nodeMap = utils::parse(newLine,',','=','"');

    // If the node is rootnode, we will create a new RootNode object and assign a reference to it in forestMap
    if ( nodeMap["NODE"] == "*" ) {
      rootNodes_[treeIdx] = new RootNode();
      forestMap[treeIdx]["*"] = rootNodes_[treeIdx];
    }
      
    Node* nodep = forestMap[treeIdx][nodeMap["NODE"]];

    // Set train prediction of the node
    nodep->setTrainPrediction(datadefs::str2num(nodeMap["PRED"]));
    
    // If the node has a splitter...
    if ( nodeMap.find("SPLITTER") != nodeMap.end() ) {
      
      size_t featureIdx = treeData_->getFeatureIdx(nodeMap["SPLITTER"]);

      // If the splitter is numerical...
      if ( nodeMap["SPLITTERTYPE"] == "NUMERICAL" ) {
       
	nodep->setSplitter(featureIdx,nodeMap["SPLITTER"],datadefs::str2num(nodeMap["LVALUES"]));
	
      } else {

	set<string> splitLeftValues = utils::keys(nodeMap["LVALUES"],':');
	set<string> splitRightValues = utils::keys(nodeMap["RVALUES"],':');

	nodep->setSplitter(featureIdx,nodeMap["SPLITTER"],splitLeftValues,splitRightValues);
      
      }

      forestMap[treeIdx][nodeMap["NODE"]+"L"] = nodep->leftChild();
      forestMap[treeIdx][nodeMap["NODE"]+"R"] = nodep->rightChild();

    }
    
  }
  
}

StochasticForest::~StochasticForest() {
  for(size_t treeIdx = 0; treeIdx < rootNodes_.size(); ++treeIdx) {
    delete rootNodes_[treeIdx];
  }

}

/* Prints the forest into a file, so that the forest can be loaded for later use (e.g. prediction).
 * This function still lacks most of the content. 
 */
void StochasticForest::printToFile(const string& fileName) {

  // Open stream for writing
  ofstream toFile( fileName.c_str() );

  if ( parameters_.model == GBT ) {
    toFile << "FOREST=GBT";
  } else if ( parameters_.model == RF ) {
    toFile << "FOREST=RF"; 
  }

  size_t targetIdx = treeData_->getFeatureIdx( targetName_ );

  toFile << ",NTREES=" << parameters_.nTrees;
  toFile << ",TARGET=" << "\"" << targetName_ << "\"";
  //toFile << ",NCATEGORIES=" << targetSupport_.size();
  toFile << ",CATEGORIES=" << "\"" << utils::join(treeData_->categories(targetIdx),',') << "\"";
  toFile << ",SHRINKAGE=" << parameters_.shrinkage << endl;

  // Save each tree in the forest
  for ( size_t treeIdx = 0; treeIdx < rootNodes_.size(); ++treeIdx ) {
    toFile << "TREE=" << treeIdx << endl;
    rootNodes_[treeIdx]->print(toFile);
  }
  
  // Close stream
  toFile.close();
  
}

void StochasticForest::learnRF() {

  // We force the shrinkage to equal to one over the number of trees
  parameters_.shrinkage = 1.0 / parameters_.nTrees;

  // If we use contrasts, we need to permute them before using
  if(parameters_.useContrasts) {
    treeData_->permuteContrasts();
  }

  //learnedModel_ = RF;

  //These parameters, and those specified in the Random Forest initiatialization, define the type of the forest generated (an RF) 
  bool sampleWithReplacement = true;
  num_t sampleSizeFraction = 1.0;
  size_t maxNodesToStop = 2 * parameters_.nMaxLeaves - 1;
  size_t minNodeSizeToStop = parameters_.nodeSize;
  bool isRandomSplit = parameters_.isRandomSplit;
  size_t nFeaturesForSplit = parameters_.mTry; //static_cast<size_t>( mTryFraction*static_cast<num_t>( treeData_->nFeatures() ) );
  bool useContrasts = parameters_.useContrasts;

  if ( nFeaturesForSplit == 0 ) {
    cerr << "StochasticForest::learnRF() -- must have at least 1 feature sampled per split" << endl;
    exit(1);
  }

  size_t targetIdx = treeData_->getFeatureIdx( targetName_ );
  size_t nCategories = targetSupport_.size();

  //Allocates memory for the root nodes
  for(size_t treeIdx = 0; treeIdx < parameters_.nTrees; ++treeIdx) {
    rootNodes_[treeIdx] = new RootNode(sampleWithReplacement,
                                       sampleSizeFraction,
                                       maxNodesToStop,
                                       minNodeSizeToStop,
                                       isRandomSplit,
                                       nFeaturesForSplit,
                                       useContrasts,
                                       nCategories);

    size_t nNodes = 1;

    RootNode::PredictionFunctionType predictionFunctionType;

    if(treeData_->isFeatureNumerical(targetIdx)) {
      predictionFunctionType = RootNode::MEAN;
    } else {
      predictionFunctionType = RootNode::MODE;
    }

    rootNodes_[treeIdx]->growTree(treeData_,
				  targetIdx,
				  predictionFunctionType,
				  oobMatrix_[treeIdx],
				  featuresInForest_[treeIdx],
				  nNodes);
    
  }

  this->updateImportanceValues();

}

void StochasticForest::learnGBT() {
  
  //parameters_.shrinkage = shrinkage;

  size_t nCategories = targetSupport_.size();

  if (nCategories > 0) {
    parameters_.nTrees *= nCategories;
  }

  // This is required since the number of trees became multiplied by the number of classes
  oobMatrix_.resize( parameters_.nTrees );
  rootNodes_.resize( parameters_.nTrees );

  //learnedModel_ = GBT;
  
  bool sampleWithReplacement = parameters_.sampleWithReplacement;
  num_t inBoxFraction = parameters_.inBoxFraction;
  size_t maxNodesToStop = 2 * parameters_.nMaxLeaves - 1;
  size_t minNodeSizeToStop = parameters_.nodeSize;
  bool isRandomSplit = parameters_.isRandomSplit;
  size_t nFeaturesForSplit = parameters_.mTry;
  bool useContrasts = parameters_.useContrasts;

  //Allocates memory for the root nodes. With all these parameters, the RootNode is now able to take full control of the splitting process
  rootNodes_.resize( parameters_.nTrees );
  for(size_t treeIdx = 0; treeIdx < parameters_.nTrees; ++treeIdx) {
    rootNodes_[treeIdx] = new RootNode(sampleWithReplacement,
                                       inBoxFraction,
                                       maxNodesToStop,
                                       minNodeSizeToStop,
                                       isRandomSplit,
                                       nFeaturesForSplit,
                                       useContrasts,
                                       nCategories);
  }
    
  if ( nCategories == 0 ) {
    this->growNumericalGBT();
  } else {
    this->growCategoricalGBT();
  }

  this->updateImportanceValues();

}

// Grow a GBT "forest" for a numerical target variable
void StochasticForest::growNumericalGBT() {

  size_t targetIdx = treeData_->getFeatureIdx( targetName_ );
  
  //A function pointer to a function "mean()" that is used to compute the node predictions with
  RootNode::PredictionFunctionType predictionFunctionType = RootNode::MEAN;

  size_t nSamples = treeData_->nSamples();
  // save a copy of the target column because it will be overwritten
  vector<num_t> trueTargetData = treeData_->getFeatureData(targetIdx);

  // Target for GBT is different for each tree
  // reference to the target column, will overwrite it
  vector<num_t> curTargetData = trueTargetData; 

  // Set the initial prediction to zero.
  vector<num_t> prediction(nSamples, 0.0);

  size_t nNodes;

  for ( size_t treeIdx = 0; treeIdx < parameters_.nTrees; ++treeIdx ) {
    // current target is the negative gradient of the loss function
    // for square loss, it is target minus current prediction
    for (size_t i = 0; i < nSamples; i++ ) {
      curTargetData[i] = trueTargetData[i] - prediction[i];
    }

    treeData_->replaceFeatureData(targetIdx,curTargetData);

    // Grow a tree to predict the current target
    rootNodes_[treeIdx]->growTree(treeData_,
				  targetIdx,
				  predictionFunctionType,
				  oobMatrix_[treeIdx],
				  featuresInForest_[treeIdx],
				  nNodes);


    // What kind of a prediction does the new tree produce?
    vector<num_t> curPrediction = StochasticForest::predictDatasetByTree(treeIdx);

    // Calculate the current total prediction adding the newly generated tree
    num_t sqErrorSum = 0.0;
    for (size_t i=0; i<nSamples; i++) {
      prediction[i] = prediction[i] + parameters_.shrinkage * curPrediction[i];

      // diagnostics
      num_t iError = trueTargetData[i]-prediction[i];
      sqErrorSum += iError*iError;

    }

  }

  // GBT-forest is now done!
  // restore true target
  treeData_->replaceFeatureData(targetIdx,trueTargetData);
}

// Grow a GBT "forest" for a categorical target variable
void StochasticForest::growCategoricalGBT() {
  
  size_t targetIdx = treeData_->getFeatureIdx( targetName_ );
  size_t nCategories = treeData_->nCategories( targetIdx );

  //A function pointer to a function "gamma()" that is used to compute the node predictions with
  RootNode::PredictionFunctionType predictionFunctionType = RootNode::GAMMA;

  // Save a copy of the target column because it will be overwritten later.
  // We also know that it must be categorical.
  size_t nSamples = treeData_->nSamples();
  vector<num_t> trueTargetData = treeData_->getFeatureData(targetIdx);
  vector<string> trueRawTargetData = treeData_->getRawFeatureData(targetIdx);
  
  // Target for GBT is different for each tree.
  // We use the original target column to save each temporary target.
  // Reference to the target column, will overwrite it:
  vector<num_t> curTargetData = trueTargetData;

  // For categorical target variables, we predict each category probability separately.
  // This target is numerical, thus we need to change the target variable type to
  // numerical from categorical.

  // Initialize class probability estimates and the predictions.
  // Note that dimensions in these two are reversed!
  vector< vector<num_t> > prediction(    nSamples, vector<num_t>( nCategories, 0.0 ) );
  vector< vector<num_t> > curPrediction( nCategories, vector<num_t>( nSamples, 0.0 ) );

  vector< vector<num_t> > curProbability( nSamples, vector<num_t>( nCategories ) );

  // Each iteration consists of numClasses_ trees,
  // each of those predicting the probability residual for each class.
  size_t numIterations = parameters_.nTrees / nCategories;

  for(size_t m=0; m < numIterations; ++m) {
    // Multiclass logistic transform of class probabilities from current probability estimates.
    for (size_t i=0; i<nSamples; i++) {
      StochasticForest::transformLogistic( prediction[i],  curProbability[i]);
      // each prediction[i] is a vector<num_t>(numClasses_)
    }

    // construct a tree for each class
    for (size_t k = 0; k < nCategories; ++k) {
      // target for class k is ...
      for (size_t i=0; i<nSamples; i++) {
        // ... the difference between true target and current prediction
        curTargetData[i] = (k==trueTargetData[i]) - curProbability[i][k];
      }

      // For each tree the target data becomes the recently computed residuals
      treeData_->replaceFeatureData(targetIdx,curTargetData);

      // Grow a tree to predict the current target
      size_t treeIdx = m * nCategories + k; // tree index
      size_t nNodes;
      rootNodes_[treeIdx]->growTree(treeData_,
				    targetIdx,
				    predictionFunctionType,
				    oobMatrix_[treeIdx],
				    featuresInForest_[treeIdx],
				    nNodes);

      //cout << "Tree " << treeIdx << " ready, predicting OOB samples..." << endl;

      // What kind of a prediction does the new tree produce
      // out of the whole training data set?
      curPrediction[k] = StochasticForest::predictDatasetByTree(treeIdx);
      // Calculate the current total prediction adding the newly generated tree
      for (size_t i = 0; i < nSamples; i++) {
        prediction[i][k] = prediction[i][k] + parameters_.shrinkage * curPrediction[k][i];
      }
    }
  }
  
  // GBT-forest is now done!
  // restore the true target
  treeData_->replaceFeatureData(targetIdx,trueRawTargetData);
}

void StochasticForest::transformLogistic(vector<num_t>& prediction, vector<num_t>& probability) {

  size_t nCategories = targetSupport_.size();

  // Multiclass logistic transform of class probabilities from current probability estimates.
  assert( nCategories == prediction.size() );
  vector<num_t>& expPrediction = probability; // just using the space by a different name

  // find maximum prediction
  vector<num_t>::iterator maxPrediction = max_element( prediction.begin(), prediction.end() );
  // scale by maximum to prevent numerical errors

  num_t expSum = 0.0;
  size_t k;
  for(k=0; k < nCategories; ++k) {
    expPrediction[k] = exp( prediction[k] - *maxPrediction ); // scale by maximum
    expSum += expPrediction[k];
  }
  for(k = 0; k < nCategories; ++k) {
    probability[k] = expPrediction[k] / expSum;
  }
}

// Use a single GBT tree to produce a prediction from a single data sample of an arbitrary data set.
num_t StochasticForest::predictSampleByTree(size_t sampleIdx, size_t treeIdx) {
    
  return( this->percolateSampleIdxByTree(sampleIdx,treeIdx)->getTrainPrediction() );
  
}  

// Use a single GBT tree to produce predictions for an arbitrary data set.
vector<num_t> StochasticForest::predictDatasetByTree(size_t treeIdx) {
  
  size_t nSamples = treeData_->nSamples();

  vector<num_t> prediction(nSamples);

  // predict for all samples
  for ( size_t sampleIdx = 0; sampleIdx < nSamples; ++sampleIdx) {
    prediction[sampleIdx] = this->predictSampleByTree(sampleIdx, treeIdx);
    // cout << "Sample " << i << ", prediction " << curPrediction[i]  << endl;
  }

  return( prediction );
}

void StochasticForest::predict(vector<string>& categoryPrediction, vector<num_t>& confidence) {

  //size_t targetIdx = treeData_->getFeatureIdx( targetName_ );
  size_t nSamples = treeData_->nSamples();
  size_t nCategories = targetSupport_.size(); //treeData_->nCategories( targetIdx );

  categoryPrediction.resize( nSamples );
  confidence.resize( nSamples );

  if ( parameters_.model == GBT ) {
    
    // For classification, each "tree" is actually numClasses_ trees in a row, each predicting the probability of its own class.
    size_t numIterations = parameters_.nTrees / nCategories;

    // Vector storing the transformed probabilities for each class prediction
    vector<num_t> prediction( nCategories );

    // Vector storing true probabilities for each class prediction
    vector<num_t> probPrediction( nCategories );

    // For each sample we need to produce predictions for each class.
    for ( size_t i = 0; i < nSamples; i++ ) {
      
      for (size_t k = 0; k < nCategories; ++k) {
	
	// Initialize the prediction
	prediction[k] = 0.0;
	
	// We go through
	for(size_t m = 0; m < numIterations; ++m) {
	  
	  // Tree index
	  size_t t =  m * nCategories + k;
	  
	  // Shrinked shift towards the new prediction
	  prediction[k] = prediction[k] + parameters_.shrinkage * this->predictSampleByTree(i, t);
	  
	}
      }
      
      // ... find index of maximum prediction, this is the predicted category
      vector<num_t>::iterator maxProb = max_element( prediction.begin(), prediction.end() );
      size_t maxProbCategory = maxProb - prediction.begin();
      
      categoryPrediction[i] = targetSupport_[maxProbCategory]; 
      //treeData_->getRawFeatureData(targetIdx,maxProbCategory); //backMapping[ maxProbCategory ]; // classes are 0,1,2,...
      
      StochasticForest::transformLogistic(prediction, probPrediction); // predictions-to-probabilities
      
      vector<num_t>::iterator largestElementIt = max_element(probPrediction.begin(),probPrediction.end());
      confidence[i] = *largestElementIt;
    }

  } else if ( parameters_.model == RF ) {

    for ( size_t i = 0; i < nSamples; ++i ) {
      
      // This container stores the categories and their frequencies
      map<string,size_t> catFrequency;

      for ( size_t t = 0; t < parameters_.nTrees; ++t ) {
	
	// Extract predicted category
	string newCategory = targetSupport_[static_cast<size_t>(this->predictSampleByTree(i, t))];
	
	// If category has been seen already, increment frequency by one
	if ( catFrequency.find(newCategory) != catFrequency.end() ) {
	  ++catFrequency[newCategory];
	} else {
	  
	  // Otherwise assign frequency to 1
	  catFrequency[newCategory] = 1;

	}
      }

      // There can be less, but no more categories represented in the support than in the data
      assert( catFrequency.size() <= nCategories );

      // Initialize confidence to 0.0
      confidence[i] = 0.0;
      
      // Loop through the categories that have predictions
      for ( map<string,size_t>::const_iterator it( catFrequency.begin() ); it != catFrequency.end(); ++it ) {
	
	// If the prediction has higher prediction than the previous best, switch  
	if ( it->second > confidence[i]*parameters_.nTrees ) {
	  categoryPrediction[i] = it->first;
	  confidence[i] = static_cast<num_t>(it->second) / parameters_.nTrees;
	}
      }

    }

  } else {

    // This should never happen
    cerr << "StochasticForest::predict(string) -- no model to predict with" << endl;
    exit(1);

  }

}

void StochasticForest::predict(vector<num_t>& prediction, vector<num_t>& confidence) {
 
  size_t nSamples = treeData_->nSamples();
  prediction.resize(nSamples);
  confidence.resize(nSamples);

  for ( size_t sampleIdx = 0; sampleIdx < nSamples; ++sampleIdx ) {
    prediction[sampleIdx] = 0;
    for ( size_t treeIdx = 0; treeIdx < parameters_.nTrees; ++treeIdx) {
      prediction[sampleIdx] = prediction[sampleIdx] + parameters_.shrinkage * predictSampleByTree(sampleIdx, treeIdx);
    }
    
  }
 
}

map<Node*,vector<size_t> > StochasticForest::percolateSampleIcsByTree(const vector<size_t>& sampleIcs, const size_t treeIdx) {
  
  // This horrible construct stores information about which
  // nodes in the tree contain which training samples
  map<Node*,vector<size_t> > trainIcs;

  // Loop through all train sample indices
  for ( size_t i = 0; i < sampleIcs.size(); ++i) {

    Node* nodep( this->percolateSampleIdxByTree(sampleIcs[i],treeIdx) );

    trainIcs[nodep].push_back(sampleIcs[i]);

  }
    
  if(false) {
    cout << "Train samples percolated accordingly:" << endl;
    size_t iter = 0;
    for(map<Node*,vector<size_t> >::const_iterator it(trainIcs.begin()); it != trainIcs.end(); ++it, ++iter) {
      cout << "leaf node " << iter << " -- prediction " << it->first->getTrainPrediction() << " :"; 
      for(size_t i = 0; i < it->second.size(); ++i) {
        cout << " " << it->second[i];
      }
      cout << endl;
    }
  }

  return( trainIcs );

}

map<Node*,vector<size_t> > StochasticForest::percolateSampleIcsByTreeAtRandom(const size_t featureIdx, const vector<size_t>& sampleIcs, const size_t treeIdx) {

  // This horrible construct stores information about which 
  // nodes in the tree contain which training samples
  map<Node*,vector<size_t> > trainIcs;

  // Loop through all train sample indices
  for ( size_t i = 0; i < sampleIcs.size(); ++i) {

    // Initialize a pointer to the root node
    Node* nodep( this->percolateSampleIdxByTreeAtRandom(featureIdx,sampleIcs[i],treeIdx) );
    
    trainIcs[nodep].push_back(sampleIcs[i]);
    
  }

  return( trainIcs );

}

Node* StochasticForest::percolateSampleIdxByTree(const size_t sampleIdx, const size_t treeIdx) {

  Node* nodep( static_cast<Node*>(rootNodes_[treeIdx]) );

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
   
    // However, if percolation could not be done 
    // (previously unobserved category or missing value), it's time to exit the loop
    if ( childNode == nodep ) {
      break;
    }

    // Update the pointer and continue percolating
    nodep = childNode;

  }

  return( nodep );

}

Node* StochasticForest::percolateSampleIdxByTreeAtRandom(const size_t featureIdx, const size_t sampleIdx, const size_t treeIdx) {
  
  Node* nodep( static_cast<Node*>(rootNodes_[treeIdx]) );

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

    if ( childNode == nodep ) {
      break;
    }

    nodep = childNode;

  }

  return( nodep );

}

// In growForest a bootstrapper was utilized to generate in-box (IB) and out-of-box (OOB) samples.
// IB samples were used to grow the forest, excluding OOB samples. In this function, these 
// previously excluded OOB samples are used for testing which features in the trained trees seem
// to contribute the most to the quality of the predictions. This is a three-fold process:
// 
// 0. for feature_i in features:
// 1. Take the original forest, percolate OOB samples across the trees all the way to the leaf nodes
//    and check how concordant the OOB and IB samples in the leafs, on average, are.
// 2. Same as with #1, but if feature_i is to make the split, sample a random value for feature_i
//    to make the split with.
// 3. Quantitate relative increase of disagreement between OOB and IB data on the leaves, with and 
//    without random sampling. Rationale: if feature_i is important, random sampling will have a 
//    big impact, thus, relative increase of disagreement will be high.  
void StochasticForest::updateImportanceValues() {

  // The number of real features in the data matrix...
  size_t nRealFeatures = treeData_->nFeatures();

  // But as there is an equal amount of contrast features, the total feature count is double that.
  size_t nAllFeatures = 2 * nRealFeatures;

  // Initialize output struct, storing the importance scores and OOB errors
  //StochasticForest::Report report;
  importanceValues_.resize( nAllFeatures, 0.0 );
  oobError_ = 0.0;

  vector<bool> hasImportance( nAllFeatures, false );
  size_t nOobSamples = 0; // !! Potentially Unintentional Humor: "Noob". That is actually intentional. :)
  
  //num_t meanTreePredictionError = 0.0;

  // The random forest object stores the mapping from trees to features it contains, which makes
  // the subsequent computations void of unnecessary looping
  for ( map<size_t, set<size_t> >::const_iterator tit(featuresInForest_.begin()); tit != featuresInForest_.end(); ++tit ) {
    size_t treeIdx = tit->first;
      
    size_t nNewOobSamples = oobMatrix_[treeIdx].size(); 
    nOobSamples += nNewOobSamples;

    map<Node*,vector<size_t> > trainIcs = this->percolateSampleIcsByTree(oobMatrix_[treeIdx],treeIdx);

    num_t treePredictionError = this->predictionError(trainIcs);

    // Accumulate tree impurity
    oobError_ += treePredictionError / parameters_.nTrees;

    assert( !datadefs::isNAN(oobError_) );

    // Loop through all features in the tree
    for ( set<size_t>::const_iterator fit( tit->second.begin()); fit != tit->second.end(); ++fit ) {
      size_t featureIdx = *fit;
    
      trainIcs = this->percolateSampleIcsByTreeAtRandom(featureIdx,oobMatrix_[treeIdx],treeIdx);
      num_t permutedTreePredictionError = this->predictionError(trainIcs);

      importanceValues_[featureIdx] += nNewOobSamples * (permutedTreePredictionError - treePredictionError);
      hasImportance[featureIdx] = true;
    }
      
  }
  
  for ( size_t featureIdx = 0; featureIdx < nAllFeatures; ++featureIdx ) {
    
    if ( !hasImportance[featureIdx] || fabs(oobError_) < datadefs::EPS ) {
      importanceValues_[featureIdx] = datadefs::NUM_NAN;
    } else {
      importanceValues_[featureIdx] /= ( nOobSamples * oobError_ );
    }
  }

  contrastImportanceValues_.resize(nRealFeatures);

  copy(importanceValues_.begin() + nRealFeatures,
       importanceValues_.end(),
       contrastImportanceValues_.begin());
  
  importanceValues_.resize(nRealFeatures);

}

/**
   Returns a vector of node counts in the trees of the forest
*/
vector<size_t> StochasticForest::nNodes() {
  
  // Get the number of trees
  size_t nTrees = this->nTrees();

  // Initialize the node count vector
  vector<size_t> nNodes(nTrees) ;

  // Loop through all trees
  for(size_t treeIdx = 0; treeIdx < nTrees; ++treeIdx) {

    // Get the node count for the tree
    nNodes[treeIdx] = this->nNodes(treeIdx);

  }
  
  return( nNodes );
}

/**
   Returns the number of nodes in tree treeIdx
*/
size_t StochasticForest::nNodes(const size_t treeIdx) {
  return( rootNodes_[treeIdx]->nNodes() );
}

map<size_t,map<size_t,size_t> > StochasticForest::featureFrequency() {

  // The output variable that stores all the frequencies 
  // of the features in the trees
  map<size_t,map<size_t,size_t> > frequency;

  // We loop through all the trees
  for ( size_t treeIdx = 0; treeIdx < parameters_.nTrees; ++treeIdx ) {
    
    for ( set<size_t>::const_iterator it1(featuresInForest_[treeIdx].begin() ); it1 != featuresInForest_[treeIdx].end(); ++it1 ) {
      set<size_t>::const_iterator it2( it1 );
      ++it2;
      while ( it2 != featuresInForest_[treeIdx].end() ) {
	
	size_t lokey,hikey;
	if ( *it1 <= *it2 ) { 
	  lokey = *it1;
	  hikey = *it2;
	} else {
	  lokey = *it2;
	  hikey = *it1;
	}

	if ( frequency.find(lokey) == frequency.end() ) {
	  frequency[lokey][hikey] = 1;
	} else {
	  ++frequency[lokey][hikey];
	}

	++it2;

      }
    }

  }
  

  return(frequency);
}

/**
   Returns the number of trees in the forest
*/
size_t StochasticForest::nTrees() {
  return( rootNodes_.size() );
}

num_t StochasticForest::predictionError(const map<Node*,vector<size_t> >& trainIcs) {

  // Get index point to the target in treeData_
  size_t targetIdx = treeData_->getFeatureIdx( targetName_ );

  // Initialize oobError to 0.0
  num_t predictionError = 0.0;

  // Count total number of samples
  size_t nSamples = 0;

  // Get the type of the target
  bool isTargetNumerical = treeData_->isFeatureNumerical(targetIdx);

  // Loop through all nodes that have samples assigned to them
  for(map<Node*,vector<size_t> >::const_iterator it(trainIcs.begin()); it != trainIcs.end(); ++it) {

    // Get target data
    // NOTE1: it->second points to the data indices
    // NOTE2: we happen to know that the indices do not point to data with NANs, which is why
    //        we don't have to make explicit NAN-checks either
    vector<num_t> targetData = treeData_->getFeatureData(targetIdx,it->second);
    
    // Get node prediction
    // NOTE: it->first points to the node
    num_t nodePrediction = it->first->getTrainPrediction();
    
    assert( !datadefs::isNAN(nodePrediction) );

    // Number of percolateds sample in node
    size_t nSamplesInNode = targetData.size();

    // Depending on the type of the target, different error calculation equation is used
    if(isTargetNumerical) {

      // Use squared error formula
      for(size_t i = 0; i < nSamplesInNode; ++i) {
        predictionError += pow(nodePrediction - targetData[i],2);
      }

    } else {
      
      // Use fraction of misprediction formula
      for(size_t i = 0; i < nSamplesInNode; ++i) {
        predictionError += nodePrediction != targetData[i]; 
      }
    }

    assert( !datadefs::isNAN(predictionError) );

    // Accumulate total sample counter
    nSamples += nSamplesInNode;

  }

  // Get mean
  if(nSamples > 0) {
    predictionError /= nSamples;
  } else {
    predictionError = datadefs::NUM_NAN;
  }

  return( predictionError );

}

vector<num_t> StochasticForest::importanceValues() {
  return( importanceValues_ );
}

vector<num_t> StochasticForest::contrastImportanceValues() {
  return( contrastImportanceValues_ );
}

num_t StochasticForest::oobError() {
  return( oobError_ );
}

