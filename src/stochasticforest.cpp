#include <cmath>
#include <ctime>
#include <iostream>
#include <iomanip>

#include "stochasticforest.hpp"
#include "datadefs.hpp"

StochasticForest::StochasticForest(Treedata* treeData, const size_t targetIdx, const size_t nTrees):
  treeData_(treeData),
  targetIdx_(targetIdx),
  nTrees_(nTrees),
  rootNodes_(0),
  oobMatrix_(0) {

  featuresInForest_.clear();

  treeData_ = treeData;
  targetIdx_ = targetIdx;

  learnedModel_ = NO_MODEL;

  //cout << "ctor: " << learnedModel_ << endl;

}

StochasticForest::~StochasticForest() {
  for(size_t treeIdx = 0; treeIdx < rootNodes_.size(); ++treeIdx) {
    delete rootNodes_[treeIdx];
  }
}

void StochasticForest::learnRF(const size_t mTry, 
			       const size_t nMaxLeaves,
			       const size_t nodeSize,
			       const bool useContrasts,
			       const bool isOptimizedNodeSplit) {

  if(useContrasts) {
    treeData_->permuteContrasts();
  }

  numClasses_ = treeData_->nCategories( targetIdx_ );
  learnedModel_ = RF_MODEL;
  featuresInForest_.clear();
  oobMatrix_.resize( nTrees_ );
  rootNodes_.resize( nTrees_ );

  //These parameters, and those specified in the Random Forest initiatialization, define the type of the forest generated (an RF) 
  bool sampleWithReplacement = true;
  num_t sampleSizeFraction = 1.0;
  size_t maxNodesToStop = 2 * nMaxLeaves - 1;
  size_t minNodeSizeToStop = nodeSize;
  bool isRandomSplit = true;
  size_t nFeaturesForSplit = mTry;

  //Allocates memory for the root nodes
  for(size_t treeIdx = 0; treeIdx < nTrees_; ++treeIdx) {
    rootNodes_[treeIdx] = new RootNode(sampleWithReplacement,
                                       sampleSizeFraction,
                                       maxNodesToStop,
                                       minNodeSizeToStop,
                                       isRandomSplit,
                                       nFeaturesForSplit,
                                       useContrasts,
                                       isOptimizedNodeSplit,
                                       numClasses_);

    size_t nNodes;

    //void (*leafPredictionFunction)(const vector<num_t>&, const size_t);
    RootNode::LeafPredictionFunctionType leafPredictionFunctionType;

    if(treeData_->isFeatureNumerical(targetIdx_)) {
      leafPredictionFunctionType = RootNode::LEAF_MEAN;
    } else {
      leafPredictionFunctionType = RootNode::LEAF_MODE;
    }

    rootNodes_[treeIdx]->growTree(treeData_,
				  targetIdx_,
				  leafPredictionFunctionType,
				  oobMatrix_[treeIdx],
				  featuresInForest_[treeIdx],
				  nNodes);
    
  }
}

void StochasticForest::learnGBT(const size_t nMaxLeaves, 
				const num_t shrinkage, 
				const num_t subSampleSize) {
  
  
  shrinkage_ = shrinkage;
  numClasses_ = treeData_->nCategories(targetIdx_);
  if (numClasses_ > 0) {
    nTrees_ *= numClasses_;
  }

  learnedModel_ = GBT_MODEL;
  featuresInForest_.clear();
  oobMatrix_.resize(nTrees_);
  rootNodes_.resize(nTrees_);
  
  bool sampleWithReplacement = false;
  num_t sampleSizeFraction = subSampleSize;
  size_t maxNodesToStop = 2 * nMaxLeaves - 1;
  size_t minNodeSizeToStop = 1;
  bool isRandomSplit = false;
  size_t nFeaturesForSplit = treeData_->nFeatures();
  bool useContrasts = false;
  bool isOptimizedNodeSplit = true; // WILL BE AN OPTION, AT THE MOMENT REGULAR SPLITTING ISN'T WORKING
  //bool isSaveLeafTrainPrediction = false;

  //Allocates memory for the root nodes. With all these parameters, the RootNode is now able to take full control of the splitting process
  rootNodes_.resize(nTrees_);
  for(size_t treeIdx = 0; treeIdx < nTrees_; ++treeIdx) {
    rootNodes_[treeIdx] = new RootNode(sampleWithReplacement,
                                       sampleSizeFraction,
                                       maxNodesToStop,
                                       minNodeSizeToStop,
                                       isRandomSplit,
                                       nFeaturesForSplit,
                                       useContrasts,
                                       isOptimizedNodeSplit,
                                       numClasses_);
  }
  //Let's grow the forest
  //cout << "Target is "<< treeData_->getFeatureName(targetIdx_) <<" ["<<targetIdx_<<"]. It has "<<numClasses_<<" classes."<<endl;
  
  if(numClasses_ == 0) {
    StochasticForest::growNumericalGBT();
  } else {
    StochasticForest::growCategoricalGBT();
  }
}

// Grow a GBT "forest" for a numerical target variable
void StochasticForest::growNumericalGBT() {
  
  //A function pointer to a function "mean()" that is used to compute the node predictions with
  //void (*leafPredictionFunction)(const vector<num_t>&, const size_t) = Node::leafMean;
  RootNode::LeafPredictionFunctionType leafPredictionFunctionType = RootNode::LEAF_MEAN;

  size_t nSamples = treeData_->nSamples();
  // save a copy of the target column because it will be overwritten
  vector<num_t> trueTargetData;
  treeData_->getFeatureData(targetIdx_,trueTargetData);

  // Target for GBT is different for each tree
  // reference to the target column, will overwrite it
  vector<num_t>& curTargetData = treeData_->features_[targetIdx_].data; //// THIS WILL BECOME INVALID UPON REMOVAL OF FRIENDSHIP ASSIGNMENT IN TREEDATA ////

  // Set the initial prediction to zero.
  vector<num_t> prediction(nSamples, 0.0);
  vector<num_t> curPrediction(nSamples);

  size_t nNodes;


  for(size_t t = 0; t < nTrees_; ++t) {
    // current target is the negative gradient of the loss function
    // for square loss, it is target minus current prediction
    for (size_t i=0; i<nSamples; i++) {
      curTargetData[i] = trueTargetData[i] - prediction[i];
    }

    // Grow a tree to predict the current target
    rootNodes_[t]->growTree(treeData_,
                            targetIdx_,
                            leafPredictionFunctionType,
                            oobMatrix_[t],
                            featuresInForest_[t],
                            nNodes);


    // What kind of a prediction does the new tree produce?
    StochasticForest::predictDatasetByTree(treeData_, t, curPrediction);

    // Calculate the current total prediction adding the newly generated tree
    //cout <<endl<<"Tree "<<t<<": Predictions:"<<endl;
    num_t sqErrorSum = 0.0;
    for (size_t i=0; i<nSamples; i++) {
      prediction[i] = prediction[i] + shrinkage_ * curPrediction[i];

      // diagnostics
      num_t iError = trueTargetData[i]-prediction[i];
      sqErrorSum += iError*iError;
      //cout << setiosflags(ios::fixed) << setprecision(4);
      ///cout <<"i="<<setw(4)<<i;
      // cout <<" cp="<<setw(8)<<curPrediction[i]<<" ct="<<setw(8)<< curTargetData[i]<<" ce="<<setw(8)<< curTargetData[i]-curPrediction[i];
      //cout << " p="<<setw(8)<<   prediction[i]<< " t="<<setw(8)<<trueTargetData[i]<< " e="<<setw(8)<< iError;
      //cout << endl;
    }
    //cout << "rmserror="<<sqrt(sqErrorSum/nSamples)<<endl;
  }
  // GBT-forest is now done!

  // restore true target column
  treeData_->features_[targetIdx_].data = trueTargetData; //// THIS WILL BECOME INVALID UPON REMOVAL OF FRIENDSHIP ASSIGNMENT IN TREEDATA ////
}

// Grow a GBT "forest" for a categorical target variable
void StochasticForest::growCategoricalGBT() {
  
  //A function pointer to a function "gamma()" that is used to compute the node predictions with
  //void (*leafPredictionFunction)(const vector<num_t>&, const size_t) = Node::leafGamma;
  RootNode::LeafPredictionFunctionType leafPredictionFunctionType = RootNode::LEAF_GAMMA;

  // Save a copy of the target column because it will be overwritten later.
  // We also know that it must be categorical.
  size_t nSamples = treeData_->nSamples();
  vector<num_t> trueTargetData;
  treeData_->getFeatureData(targetIdx_,trueTargetData);

  // Target for GBT is different for each tree.
  // We use the original target column to save each temporary target.
  // Reference to the target column, will overwrite it:
  vector<num_t>& curTargetData = treeData_->features_[targetIdx_].data; //// THIS WILL BECOME INVALID UPON REMOVAL OF FRIENDSHIP ASSIGNMENT IN TREEDATA ////

  // For categorical target variables, we predict each category probability separately.
  // This target is numerical, thus we need to change the target variable type to
  // numerical from categorical.
  treeData_->features_[targetIdx_].isNumerical = true; //// THIS WILL BECOME INVALID UPON REMOVAL OF FRIENDSHIP ASSIGNMENT IN TREEDATA ////


  // Initialize class probability estimates and the predictions.
  // Note that dimensions in these two are reversed!
  vector< vector<num_t> > prediction(    nSamples, vector<num_t>( numClasses_, 0.0 ) );
  vector< vector<num_t> > curPrediction( numClasses_, vector<num_t>( nSamples, 0.0 ) );

  vector< vector<num_t> > curProbability( nSamples, vector<num_t>( numClasses_)  );
  // TODO don't really need all this space allocated
  //vector<size_t> oobIcs;
  //set<size_t> featuresInTree;
  //bool useContrasts = false;

  // Each iteration consists of numClasses_ trees,
  // each of those predicting the probability residual for each class.
  size_t numIterations = nTrees_ / numClasses_;

  for(size_t m=0; m < numIterations; ++m) {
    // Multiclass logistic transform of class probabilities from current probability estimates.
    for (size_t i=0; i<nSamples; i++) {
      StochasticForest::transformLogistic( numClasses_, prediction[i],  curProbability[i]);
      // each prediction[i] is a vector<num_t>(numClasses_)
    }

    // construct a tree for each class
    for (size_t k = 0; k < numClasses_; ++k) {
      // target for class k is ...
      for (size_t i=0; i<nSamples; i++) {
        // ... the difference between true target and current prediction
        curTargetData[i] = (k==trueTargetData[i]) - curProbability[i][k];
      }

      // Grow a tree to predict the current target
      size_t t = m * numClasses_ + k; // tree index
      size_t nNodes;
      rootNodes_[t]->growTree(treeData_,
                              targetIdx_,
                              leafPredictionFunctionType,
                              oobMatrix_[t],
                              featuresInForest_[t],
                              nNodes);

      // What kind of a prediction does the new tree produce
      // out of the whole training data set?
      StochasticForest::predictDatasetByTree(treeData_, t, curPrediction[k] );
      // Calculate the current total prediction adding the newly generated tree
      //cout <<"Iter="<<m<<" Class="<<k<<": Predictions:"<<endl;
      for (size_t i=0; i<nSamples; i++) {
        prediction[i][k] = prediction[i][k] + shrinkage_ * curPrediction[k][i];
        // cout <<"i="<<i<<" CurPrediction="<<curPrediction[k][i]<<" prediction="<<prediction[i][k]<<endl;
      }
    }
  }

  // GBT-forest is now done!
  // restore the true target column and its true type
  treeData_->features_[targetIdx_].data = trueTargetData; //// THIS WILL BECOME INVALID UPON REMOVAL OF FRIENDSHIP ASSIGNMENT IN TREEDATA ////
  treeData_->features_[targetIdx_].isNumerical = false; //// THIS WILL BECOME INVALID UPON REMOVAL OF FRIENDSHIP ASSIGNMENT IN TREEDATA ////
}

void StochasticForest::transformLogistic(const size_t numClasses, vector<num_t>& prediction, vector<num_t>& probability) {
  // Multiclass logistic transform of class probabilities from current probability estimates.
  assert( numClasses_ == prediction.size() );
  vector<num_t>& expPrediction = probability; // just using the space by a different name

  // find maximum prediction
  vector<num_t>::iterator maxPrediction = max_element( prediction.begin(), prediction.end() );
  // scale by maximum to prevent numerical errors

  num_t expSum = 0.0;
  size_t k;
  for(k=0; k < numClasses_; ++k) {
    expPrediction[k] = exp( prediction[k] - *maxPrediction ); // scale by maximum
    expSum += expPrediction[k];
  }
  for(k = 0; k < numClasses_; ++k) {
    probability[k] = expPrediction[k] / expSum;
  }
}

// Use a single GBT tree to produce a prediction from a single data sample of an arbitrary data set.
num_t StochasticForest::predictSampleByTree(Treedata* treeData, size_t sampleIdx, size_t treeIdx) {
  
  // Root of current tree
  Node *currentNode = rootNodes_[treeIdx];
  
  // Traverse to the leaf of the tree
  while(currentNode->hasChildren()) {
    
    num_t value;
    
    // Get the splitter of the branch point
    size_t featureIdx = currentNode->getSplitter();

    // Get the value of the splitter feature of the chosen sample 
    treeData->getFeatureData(featureIdx, sampleIdx, value);
    
    // The node then makes the branch decision. The chosen child node becomes the new currentNode 
    currentNode = currentNode->percolateData(value);
  }

  // The loop has ended, and currentNode now points to a leaf node; get the prediction
  return currentNode->getLeafTrainPrediction();
}

// Use a single GBT tree to produce predictions for an arbitrary data set.
void StochasticForest::predictDatasetByTree(Treedata* treeData, size_t treeIdx, vector<num_t>& curPrediction) {
  size_t nSamples = treeData->nSamples();
  // predict for all samples
  for(size_t i = 0; i < nSamples; ++i) {
    curPrediction[i] = StochasticForest::predictSampleByTree(treeData, i, treeIdx);
    // cout << "Sample " << i << ", prediction " << curPrediction[i]  << endl;
  }
}

void StochasticForest::predict(vector<string>& prediction, vector<num_t>& confidence) {
  assert(numClasses_ != 0);
  StochasticForest::predict(treeData_, prediction, confidence);
}

void StochasticForest::predict(vector<num_t>& prediction, vector<num_t>& confidence) {
  assert(numClasses_ == 0);
  StochasticForest::predict(treeData_, prediction, confidence);
}

// Predict with the trained model and using an arbitrary data set. StochasticForest::numClasses_ determines
// whether the prediction is for a categorical or a numerical variable.
void StochasticForest::predict(Treedata* treeData, vector<string>& prediction, vector<num_t>& confidence) {
  assert( numClasses_ != 0 );

  switch ( learnedModel_ ) {
  case GBT_MODEL:
    StochasticForest::predictWithCategoricalGBT(treeData, prediction, confidence);
    break;
  case RF_MODEL:
    cerr << "Implementation of categorical prediction with RFs isn't yet working" << endl;
    assert(false);
    StochasticForest::predictWithCategoricalRF(treeData, prediction);
    break;
  case NO_MODEL:
    cerr << "Cannot predict -- no model trained" << endl;
    assert(false);
  }

}

// Predict with the trained model and using an arbitrary data set. StochasticForest::numClasses_ determines
// whether the prediction is for a categorical or a numerical variable.
void StochasticForest::predict(Treedata* treeData, vector<num_t>& prediction, vector<num_t>& confidence) {
  assert( numClasses_ == 0 );

  //cout << "Prediction of numerical data is unstable and doesn't yet yield confidence metrics" << endl;
  
  switch ( learnedModel_ ) {
  case GBT_MODEL:
    StochasticForest::predictWithNumericalGBT(treeData, prediction, confidence);
    break;
  case RF_MODEL:
    StochasticForest::predictWithNumericalRF(treeData, prediction);
    break;
  case NO_MODEL:
    cerr << "Cannot predict -- no model trained" << endl;
    assert(false);
  }

}

void StochasticForest::predictWithCategoricalRF(Treedata* treeData, vector<string>& categoryPrediction) {

  cerr << "Prediction with RF isn't yet working" << endl;
  assert(false);

  categoryPrediction.resize( treeData->nSamples() );
  

}

void StochasticForest::predictWithNumericalRF(Treedata* treeData, vector<num_t>& prediction) {

  cerr << "Prediction with RF isn't yet working" << endl;
  assert(false);

  prediction.resize( treeData->nSamples() );

  for(size_t sampleIdx = 0; sampleIdx < treeData->nSamples(); ++sampleIdx) {
    
    prediction[sampleIdx] = 0.0;
    
    for(size_t treeIdx = 0; treeIdx < nTrees_; ++treeIdx) {
      
      Node* nodep(rootNodes_[treeIdx]);
      StochasticForest::percolateSampleIdx(treeData, sampleIdx, &nodep);
      prediction[sampleIdx] += nodep->getLeafTrainPrediction();
    
    }

    prediction[sampleIdx] /= nTrees_;
  
  }
  
}


// Predict categorical target using a GBT "forest" from an arbitrary data set. The function also outputs a confidence score for 
// the predictions. In this case the confidence score is the probability for the prediction.
void StochasticForest::predictWithCategoricalGBT(Treedata* treeData, vector<string>& categoryPrediction, vector<num_t>& confidence) {
  
  size_t nSamples = treeData->nSamples();
  
  // For classification, each "tree" is actually numClasses_ trees in a row, each predicting the probability of its own class.
  size_t numIterations = nTrees_ / numClasses_;

  // Vector storing the transformed probabilities for each class prediction
  vector<num_t> prediction( numClasses_ );

  // Vector storing true probabilities for each class prediction
  vector<num_t> probPrediction( numClasses_ );

  categoryPrediction.resize( nSamples );
  confidence.resize( nSamples );

  // For each sample we need to produce predictions for each class.
  for ( size_t i = 0; i < nSamples; i++ ) { 
    
    for (size_t k = 0; k < numClasses_; ++k) { 
      
      // Initialize the prediction
      prediction[k] = 0.0;
      
      // We go through 
      for(size_t m = 0; m < numIterations; ++m) {
      
	// Tree index
	size_t t =  m * numClasses_ + k;

	// Shrinked shift towards the new prediction
        prediction[k] = prediction[k] + shrinkage_ * predictSampleByTree(treeData, i, t);
      
      }
    }

    // ... find index of maximum prediction, this is the predicted category
    vector<num_t>::iterator maxProb = max_element( prediction.begin(), prediction.end() );
    num_t maxProbCategory = 1.0*(maxProb - prediction.begin()); 
    
    map<num_t,string> backMapping = treeData->features_[targetIdx_].backMapping;
    categoryPrediction[i] = backMapping[ maxProbCategory ]; // classes are 0,1,2,...
    
    StochasticForest::transformLogistic(numClasses_, prediction, probPrediction); // predictions-to-probabilities
    
    vector<num_t>::iterator largestElementIt = max_element(probPrediction.begin(),probPrediction.end());
    confidence[i] = *largestElementIt;
  }
}

// Predict numerical target using a GBT "forest" from an arbitrary data set
void StochasticForest::predictWithNumericalGBT(Treedata* treeData, vector<num_t>& prediction, vector<num_t>& confidence) {

  size_t nSamples = treeData->nSamples();
  prediction.resize(nSamples);
  confidence.resize(nSamples);

  //cout << "Predicting "<<nSamples<<" samples. Target="<<targetIdx_<<endl;
  for (size_t i=0; i<nSamples; i++) {
    prediction[i] = 0;
    for(size_t t = 0; t < nTrees_; ++t) {
      prediction[i] = prediction[i] + shrinkage_ * predictSampleByTree(treeData, i, t);
    }
    // diagnostic print out the true and the prediction
    //cout << i << "\t" << treeData_->features_[targetIdx_].data[i] << "\t" << prediction[i] <<endl; //// THIS WILL BECOME INVALID UPON REMOVAL OF FRIENDSHIP ASSIGNMENT IN TREEDATA ////
	}
}


void StochasticForest::percolateSampleIcs(Treedata* treeData, Node* rootNode, const vector<size_t>& sampleIcs, map<Node*,vector<size_t> >& trainIcs) {
  
  trainIcs.clear();
  //map<Node*,vector<size_t> > trainics;
  
  for(size_t i = 0; i < sampleIcs.size(); ++i) {
    Node* nodep(rootNode);
    size_t sampleIdx = sampleIcs[i];
    StochasticForest::percolateSampleIdx(treeData,sampleIdx,&nodep);
    map<Node*,vector<size_t> >::iterator it(trainIcs.find(nodep));
    if(it == trainIcs.end()) {
      Node* foop(nodep);
      vector<size_t> foo(1);
      foo[0] = sampleIdx;
      trainIcs.insert(pair<Node*,vector<size_t> >(foop,foo));
    } else {
      trainIcs[it->first].push_back(sampleIdx);
    }
      
  }
  
  
  if(false) {
    cout << "Train samples percolated accordingly:" << endl;
    size_t iter = 0;
    for(map<Node*,vector<size_t> >::const_iterator it(trainIcs.begin()); it != trainIcs.end(); ++it, ++iter) {
      cout << "leaf node " << iter << ":"; 
      for(size_t i = 0; i < it->second.size(); ++i) {
        cout << " " << it->second[i];
      }
      cout << endl;
    }
  }
}

void StochasticForest::percolateSampleIcsAtRandom(Treedata* treeData, const size_t featureIdx, Node* rootNode, const vector<size_t>& sampleIcs, map<Node*,vector<size_t> >& trainIcs) {

  trainIcs.clear();

  for(size_t i = 0; i < sampleIcs.size(); ++i) {
    Node* nodep(rootNode);
    size_t sampleIdx = sampleIcs[i];
    StochasticForest::percolateSampleIdxAtRandom(treeData,featureIdx,sampleIdx,&nodep);
    map<Node*,vector<size_t> >::iterator it(trainIcs.find(nodep));
    if(it == trainIcs.end()) {
      Node* foop(nodep);
      vector<size_t> foo(1);
      foo[0] = sampleIdx;
      trainIcs.insert(pair<Node*,vector<size_t> >(foop,foo));
    } else {
      trainIcs[it->first].push_back(sampleIdx);
    }

  }
}

void StochasticForest::percolateSampleIdx(Treedata* treeData, const size_t sampleIdx, Node** nodep) {
  while((*nodep)->hasChildren()) {
    int featureIdxNew((*nodep)->getSplitter());
    num_t value;
    treeData->getFeatureData(featureIdxNew,sampleIdx,value);
    *nodep = (*nodep)->percolateData(value);
  }
}

void StochasticForest::percolateSampleIdxAtRandom(Treedata* treeData, const size_t featureIdx, const size_t sampleIdx, Node** nodep) {
  while((*nodep)->hasChildren()) {
    size_t featureIdxNew = (*nodep)->getSplitter();
    num_t value = datadefs::NUM_NAN;
    if(featureIdx == featureIdxNew) {
      while(datadefs::isNAN(value)) {
        treeData->getRandomData(featureIdxNew,value);
      }
    } else {
      treeData->getFeatureData(featureIdxNew,sampleIdx,value);
    }
    *nodep = (*nodep)->percolateData(value);
  }
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
vector<num_t> StochasticForest::featureImportance() {

  // The number of real features in the data matrix...
  size_t nRealFeatures = treeData_->nFeatures();

  // But as there is an equal amount of contrast features, the total feature count is double that.
  size_t nAllFeatures = 2 * nRealFeatures;

  //Each feature, either real or contrast, will have a slot into which the importance value will be put.
  vector<num_t> importance( nAllFeatures, 0.0 );
  size_t nOobSamples = 0; // !! Potentially Unintentional Humor: "Noob". That is actually intentional. :)
  
  size_t nContrastsInForest = 0;

  // The random forest object stores the mapping from trees to features it contains, which makes
  // the subsequent computations void of unnecessary looping
  for ( map<size_t, set<size_t> >::const_iterator tit(featuresInForest_.begin()); tit != featuresInForest_.end(); ++tit ) {
    size_t treeIdx = tit->first;
      
    size_t nNewOobSamples = oobMatrix_[treeIdx].size(); 
    nOobSamples += nNewOobSamples;

    map<Node*,vector<size_t> > trainIcs;
    StochasticForest::percolateSampleIcs(treeData_,rootNodes_[treeIdx],oobMatrix_[treeIdx],trainIcs);
    num_t treeImpurity;
    StochasticForest::treeImpurity(treeData_,trainIcs,treeImpurity);
    //cout << "#nodes_with_train_samples=" << trainics.size() << endl;

    for ( set<size_t>::const_iterator fit( tit->second.begin()); fit != tit->second.end(); ++fit ) {
      size_t featureIdx = *fit;
    
      if ( featureIdx >= nRealFeatures ) {
        ++nContrastsInForest;
      }

      StochasticForest::percolateSampleIcsAtRandom(treeData_,featureIdx,rootNodes_[treeIdx],oobMatrix_[treeIdx],trainIcs);
      num_t permutedTreeImpurity;
      StochasticForest::treeImpurity(treeData_,trainIcs,permutedTreeImpurity);
      if ( fabs( treeImpurity ) > datadefs::EPS) {
        importance[featureIdx] += nNewOobSamples * (permutedTreeImpurity - treeImpurity) / treeImpurity;
      }
      //cout << treeData_->getFeatureName(featureIdx) << " += " << importance[featureIdx] << endl;
    }
      
  }
  
  for ( size_t featureIdx = 0; featureIdx < nAllFeatures; ++featureIdx ) {
    //importance[featureIdx] *= 100.0*nTrees_/nNodesInForest;//1.0*nNodesInForest/nTrees_; //nContrastsInForest
    //cout << "I(" << treeData_->getFeatureName(featureIdx) << ") = " << importance[featureIdx] << endl;
    if ( fabs( importance[featureIdx] ) < datadefs::EPS ) {
      importance[featureIdx] = datadefs::NUM_NAN;
    } else {
      importance[featureIdx] /= nOobSamples;
    }
    //cout << "I(" << treeData_->getFeatureName(featureIdx) << ") = " << importance[featureIdx] << endl;
  }

  return(importance);

}

size_t StochasticForest::nNodes() {
  size_t nNodes = 0;
  for(size_t treeIdx = 0; treeIdx < rootNodes_.size(); ++treeIdx) {
    nNodes += rootNodes_[treeIdx]->nNodes();
  }
  
  return(nNodes);
}

vector<size_t> StochasticForest::featureFrequency() {
  assert(false);
  vector<size_t> frequency(0);
  return(frequency);
}


void StochasticForest::treeImpurity(Treedata* treeData, 
				    map<Node*,vector<size_t> >& trainIcs, 
				    num_t& impurity) {

  impurity = 0.0;
  size_t n_tot = 0;

  bool isTargetNumerical = treeData->isFeatureNumerical(targetIdx_);

  for(map<Node*,vector<size_t> >::iterator it(trainIcs.begin()); it != trainIcs.end(); ++it) {

    vector<num_t> targetData;
    treeData->getFeatureData(targetIdx_,it->second,targetData);
    num_t leafPrediction = it->first->getLeafTrainPrediction();
    num_t leafImpurity = 0;
    size_t nSamplesInLeaf = targetData.size();
      
    if(isTargetNumerical) {
      for(size_t i = 0; i < nSamplesInLeaf; ++i) {
        leafImpurity += pow(leafPrediction - targetData[i],2);
      }

    } else {
      for(size_t i = 0; i < nSamplesInLeaf; ++i) {
        leafImpurity += leafPrediction != targetData[i]; 
	//{
	// ++leafImpurity;
        //}
      }
    }

    n_tot += nSamplesInLeaf;
    impurity += nSamplesInLeaf * leafImpurity;
  }

  if(n_tot > 0) {
    impurity /= n_tot;
  }
}
