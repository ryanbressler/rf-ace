#include "GBT.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>

GBT::GBT(Treedata* treeData, size_t targetIdx, size_t nTrees, size_t nMaxLeaves, num_t shrinkage, num_t subSampleSize):
  treeData_(treeData),
  targetIdx_(targetIdx),
  nTrees_(nTrees),
  nMaxLeaves_(nMaxLeaves),
  shrinkage_(shrinkage),
  subSampleSize_(subSampleSize) {
  
  nodeSize_ = 1;   // Allow split down to individual example
  // Check the type of the target.
  numClasses_ = treeData_->nCategories(targetIdx);
  if (numClasses_ > 0)  {
    nTrees_*= numClasses_;
  }
  
  bool sampleWithReplacement = false;
  num_t sampleSizeFraction = subSampleSize_;
  size_t maxNodesToStop = 2 * nMaxLeaves_ - 1;
  size_t minNodeSizeToStop = nodeSize_;
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
  cout << "Target is "<< treeData_->getFeatureName(targetIdx_) <<" ["<<targetIdx_<<"]. It has "<<numClasses_<<" classes."<<endl;
  GBT::growForest();
  }

GBT::~GBT() {
  for(size_t treeIdx = 0; treeIdx < nTrees_; ++treeIdx) {
    delete rootNodes_[treeIdx];
  }
}



// Grow a GBT "forest"
void GBT::growForest() {
  if ( 0 == numClasses_ )
    GBT::growForestNumerical();
  else
    GBT::growForestCategorical();
}



// Grow a GBT "forest" for a numerical target variable
void GBT::growForestNumerical() {

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

  //bool sampleWithReplacement = false;
  //num_t sampleSize = subSampleSize_;
  //size_t maxNodesToStop = 2*nMaxLeaves_ - 1;
  //size_t minNodeSizeToStop = 1;
  //bool isRandomSplit = false;
  //size_t nFeaturesInSample = treeData_->nFeatures();
  size_t nNodes;
  vector<size_t> oobIcs;
  set<size_t> featuresInTree;
  //bool useContrasts = false;

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
                            oobIcs,
                            featuresInTree,
                            nNodes);


    // What kind of a prediction does the new tree produce?
    GBT::predictDatasetByTree(t, curPrediction);

    // Calculate the current total prediction adding the newly generated tree
    cout <<endl<<"Tree "<<t<<": Predictions:"<<endl;
    num_t sqErrorSum = 0.0;
    for (size_t i=0; i<nSamples; i++) {
      prediction[i] = prediction[i] + shrinkage_ * curPrediction[i];

      // diagnostics
      num_t iError = trueTargetData[i]-prediction[i];
      sqErrorSum += iError*iError;
      cout << setiosflags(ios::fixed) << setprecision(4);
      cout <<"i="<<setw(4)<<i;
      cout <<" cp="<<setw(8)<<curPrediction[i]<<" ct="<<setw(8)<< curTargetData[i]<<" ce="<<setw(8)<< curTargetData[i]-curPrediction[i];
      cout << " p="<<setw(8)<<   prediction[i]<< " t="<<setw(8)<<trueTargetData[i]<< " e="<<setw(8)<< iError;
      cout << endl;
    }
    cout << "rmserror="<<sqrt(sqErrorSum/nSamples)<<endl;
  }
  // GBT-forest is now done!

  // restore true target column
  treeData_->features_[targetIdx_].data = trueTargetData; //// THIS WILL BECOME INVALID UPON REMOVAL OF FRIENDSHIP ASSIGNMENT IN TREEDATA ////
}



// Grow a GBT "forest" for a categorical target variable
void GBT::growForestCategorical() {

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

  //bool sampleWithReplacement = false;
  //num_t sampleSize = subSampleSize_;
  //size_t maxNodesToStop = 2*nMaxLeaves_ - 1;
  //size_t minNodeSizeToStop = 1;
  //bool isRandomSplit = false;
  //size_t nFeaturesInSample = treeData_->nFeatures();
  size_t nNodes;
  vector<size_t> oobIcs;
  set<size_t> featuresInTree;
  //bool useContrasts = false;

  // Each iteration consists of numClasses_ trees,
  // each of those predicting the probability residual for each class.
  size_t numIterations = nTrees_ / numClasses_;

  for(size_t m=0; m < numIterations; ++m) {
    // Multiclass logistic transform of class probabilities from current probability estimates.
    for (size_t i=0; i<nSamples; i++) {
      GBT::transformLogistic( prediction[i],  curProbability[i]);
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
      rootNodes_[t]->growTree(treeData_,
                              targetIdx_,
                              leafPredictionFunctionType,
                              oobIcs,
                              featuresInTree,
                              nNodes);

      // What kind of a prediction does the new tree produce
      // out of the whole training data set?
      GBT::predictDatasetByTree(t, curPrediction[k] );

      // Calculate the current total prediction adding the newly generated tree
      cout <<"Iter="<<m<<" Class="<<k<<": Predictions:"<<endl;
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



void GBT::transformLogistic(vector<num_t>& prediction, vector<num_t>& probability) {
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



// Predict using a GBT "forest" from an arbitrary data set. GBT::numClasses_ determines
// whether the prediction is for a categorical or a numerical variable.
void GBT::predictForest(Treedata* treeData, vector<num_t>& prediction) {
  if ( 0 == numClasses_ )
    predictForestNumerical(treeData, prediction);
  else
    predictForestCategorical(treeData, prediction);
}


// Predict categorical target using a GBT "forest" from an arbitrary data set.
void GBT::predictForestCategorical(Treedata* treeData, vector<num_t>& categoryPrediction) {
  size_t nSamples = treeData->nSamples();
  cout << "Predicting "<<nSamples<<" samples. Target="<<targetIdx_<<". Classes="<<numClasses_<<endl;

  // For classification, each "tree" is actually numClasses_ trees in a row
  // each predicting the probability of its own class.
  size_t numIterations = nTrees_ / numClasses_;
  size_t errorCnt = 0;
  vector<num_t> probPrediction( numClasses_ );
  vector<num_t> prediction( numClasses_ );

  for (size_t i=0; i<nSamples; i++) { // for each sample we need to ...
    for (size_t k=0; k<numClasses_; ++k) { // ... produce predictions for each class, and then ...
      prediction[k] = 0;
      for(size_t m = 0; m < numIterations; ++m) {
        size_t t =  m*numClasses_ + k; // tree index
        prediction[k] = prediction[k] + shrinkage_ * predictSampleByTree(i, t);
      }
    }

    // ... find index of maximum prediction, this is the predicted cateory
    vector<num_t>::iterator maxProb = max_element( prediction.begin(), prediction.end() );
    categoryPrediction[i] = maxProb - prediction.begin() ; // classes are 0,1,2,...

    // diagnostic print out the true and the prediction
    errorCnt += (treeData->features_[targetIdx_].data[i] != categoryPrediction[i]); //// THIS WILL BECOME INVALID UPON REMOVAL OF FRIENDSHIP ASSIGNMENT IN TREEDATA ////
    cout << i << "\t" << treeData->features_[targetIdx_].data[i] << "\t" << categoryPrediction[i]; //// THIS WILL BECOME INVALID UPON REMOVAL OF FRIENDSHIP ASSIGNMENT IN TREEDATA ////
    GBT::transformLogistic(prediction, probPrediction); // predictions-to-probabilities
    for (size_t k=0; k<numClasses_; ++k) {
      cout <<"\t"<<probPrediction[k];
    }
    cout <<endl;
  }
  cout<<"Error "<<errorCnt<<"/"<<nSamples<<"="<<100*errorCnt/nSamples<<"%"<<endl;
}



// Predict numerical target using a GBT "forest" from an arbitrary data set
void GBT::predictForestNumerical(Treedata* treeData, vector<num_t>& prediction) {
  size_t nSamples = treeData->nSamples();
  cout << "Predicting "<<nSamples<<" samples. Target="<<targetIdx_<<endl;
  for (size_t i=0; i<nSamples; i++) {
    prediction[i] = 0;
    for(size_t t = 0; t < nTrees_; ++t) {
      prediction[i] = prediction[i] + shrinkage_ * predictSampleByTree(i, t);
    }
    // diagnostic print out the true and the prediction
    cout << i << "\t" << treeData_->features_[targetIdx_].data[i] << "\t" << prediction[i] <<endl; //// THIS WILL BECOME INVALID UPON REMOVAL OF FRIENDSHIP ASSIGNMENT IN TREEDATA ////
  }
}



// Use a single GBT tree to produce a prediction from a single data sample
// of an arbitrary data set. This wouldn't have to be the train set. // IF we use another treeData_ !!
num_t GBT::predictSampleByTree(size_t sampleIdx, size_t treeIdx) {
  // root of current tree
  Node *currentNode = rootNodes_[treeIdx]; // OLD: &forest_[treeIdx][0];
  num_t value;
  // traverse to the leaf of the tree
  while(currentNode->hasChildren()) {
    size_t featureIdx = currentNode->getSplitter();
    treeData_->getFeatureData(featureIdx, sampleIdx, value);
    currentNode = currentNode->percolateData(value);
  }
  // now currentNode points to a leaf node: get the prediction
  return currentNode->getLeafTrainPrediction();
}



// Use a single GBT tree to produce predictions for a whole training data set
void GBT::predictDatasetByTree(size_t treeIdx, vector<num_t>& curPrediction) {
  size_t nSamples = treeData_->nSamples();
  // predict for all samples
  for(size_t i = 0; i < nSamples; ++i) {
    curPrediction[i] = GBT::predictSampleByTree(i, treeIdx);
    // cout << "Sample " << i << ", prediction " << curPrediction[i]  << endl;
  }
}

