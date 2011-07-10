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
  subSampleSize_(subSampleSize)
{
  nodeSize_ = 1;   // Allow split down to individual example
  // Check the type of the target.
  numClasses_ = treeData_->nCategories(targetIdx);
  if (numClasses_ > 0) 
  {
    nTrees_*= numClasses_;
  }
  
  // Allocate memory for the root nodes
  rootNodes_ = vector<Node *>( nTrees_ );
  for(size_t treeIdx = 0; treeIdx < nTrees_; ++treeIdx)
    {
      rootNodes_[treeIdx] = new Node;
    }

  //Let's grow the forest
  cout << "Target is "<< treeData_->getFeatureName(targetIdx_) <<" ["<<targetIdx_<<"]. It has "<<numClasses_<<" classes."<<endl;
  GBT::growForest();
}




GBT::~GBT()
{
  for(size_t treeIdx = 0; treeIdx < nTrees_; ++treeIdx)
    {
      delete rootNodes_[treeIdx];
    }
}



// Grow a GBT "forest"
void GBT::growForest()
{
  if ( 0 == numClasses_ )
    GBT::growForestNumerical();
  else
    GBT::growForestCategorical();
}



// Grow a GBT "forest" for a numerical target variable
void GBT::growForestNumerical()
{
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

  for(size_t t = 0; t < nTrees_; ++t)
    {
      // current target is the negative gradient of the loss function
      // for square loss, it is target minus current prediction
      for (size_t i=0; i<nSamples; i++)
        {
          curTargetData[i] = trueTargetData[i] - prediction[i];
        }

      // Grow a tree to predict the current target
      GBT::growTree(t);

      // What kind of a prediction does the new tree produce?
      GBT::predictDatasetByTree(t, curPrediction);

      // Calculate the current total prediction adding the newly generated tree
      cout <<endl<<"Tree "<<t<<": Predictions:"<<endl;
      num_t sqErrorSum = 0.0;
      for (size_t i=0; i<nSamples; i++)
        {
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
void GBT::growForestCategorical()
{
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


  // Each iteration consists of numClasses_ trees,
  // each of those predicting the probability residual for each class.
  size_t numIterations = nTrees_ / numClasses_;

  for(size_t m=0; m < numIterations; ++m)
    {
      // Multiclass logistic transform of class probabilities from current probability estimates.
      for (size_t i=0; i<nSamples; i++)
        {
          GBT::transformLogistic( prediction[i],  curProbability[i]);
          // each prediction[i] is a vector<num_t>(numClasses_)
        }

      // construct a tree for each class
      for (size_t k = 0; k < numClasses_; ++k)
        {
          // target for class k is ...
          for (size_t i=0; i<nSamples; i++)
            {
              // ... the difference between true target and current prediction
              curTargetData[i] = (k==trueTargetData[i]) - curProbability[i][k];
            }

          // Grow a tree to predict the current target
          size_t t = m * numClasses_ + k; // tree index
          GBT::growTree(t);

          // What kind of a prediction does the new tree produce
          // out of the whole training data set?
          GBT::predictDatasetByTree(t, curPrediction[k] );

          // Calculate the current total prediction adding the newly generated tree
          cout <<"Iter="<<m<<" Class="<<k<<": Predictions:"<<endl;
          for (size_t i=0; i<nSamples; i++)
            {
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



void GBT::transformLogistic(vector<num_t>& prediction, vector<num_t>& probability)
{
  // Multiclass logistic transform of class probabilities from current probability estimates.
  assert( numClasses_ == prediction.size() );
  vector<num_t>& expPrediction = probability; // just using the space by a different name

  // find maximum prediction
  vector<num_t>::iterator maxPrediction = max_element( prediction.begin(), prediction.end() );
  // scale by maximum to prevent numerical errors

  num_t expSum = 0.0;
  size_t k;
  for(k=0; k < numClasses_; ++k)
    {
      expPrediction[k] = exp( prediction[k] - *maxPrediction ); // scale by maximum
      expSum += expPrediction[k];
    }
  for(k = 0; k < numClasses_; ++k)
    {
      probability[k] = expPrediction[k] / expSum;
    }
}



// Predict using a GBT "forest" from an arbitrary data set. GBT::numClasses_ determines
// whether the prediction is for a categorical or a numerical variable.
void GBT::predictForest(Treedata* treeData, vector<num_t>& prediction)
{
  if ( 0 == numClasses_ )
    predictForestNumerical(treeData, prediction);
  else
    predictForestCategorical(treeData, prediction);
}


// Predict categorical target using a GBT "forest" from an arbitrary data set.
void GBT::predictForestCategorical(Treedata* treeData, vector<num_t>& categoryPrediction)
{
  size_t nSamples = treeData->nSamples();
  cout << "Predicting "<<nSamples<<" samples. Target="<<targetIdx_<<". Classes="<<numClasses_<<endl;

  // For classification, each "tree" is actually numClasses_ trees in a row
  // each predicting the probability of its own class.
  size_t numIterations = nTrees_ / numClasses_;
  size_t errorCnt = 0;
  vector<num_t> probPrediction( numClasses_ );
  vector<num_t> prediction( numClasses_ );

  for (size_t i=0; i<nSamples; i++) // for each sample we need to ...
    {
      for (size_t k=0; k<numClasses_; ++k) // ... produce predictions for each class, and then ...
        {
          prediction[k] = 0;
          for(size_t m = 0; m < numIterations; ++m)
            {
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
      for (size_t k=0; k<numClasses_; ++k)
        {
          cout <<"\t"<<probPrediction[k];
        }
      cout <<endl;
    }
  cout<<"Error "<<errorCnt<<"/"<<nSamples<<"="<<100*errorCnt/nSamples<<"%"<<endl;
}



// Predict numerical target using a GBT "forest" from an arbitrary data set
void GBT::predictForestNumerical(Treedata* treeData, vector<num_t>& prediction)
{
  size_t nSamples = treeData->nSamples();
  cout << "Predicting "<<nSamples<<" samples. Target="<<targetIdx_<<endl;
  for (size_t i=0; i<nSamples; i++)
    {
      prediction[i] = 0;
      for(size_t t = 0; t < nTrees_; ++t)
        {
          prediction[i] = prediction[i] + shrinkage_ * predictSampleByTree(i, t);
        }
      // diagnostic print out the true and the prediction
      cout << i << "\t" << treeData_->features_[targetIdx_].data[i] << "\t" << prediction[i] <<endl; //// THIS WILL BECOME INVALID UPON REMOVAL OF FRIENDSHIP ASSIGNMENT IN TREEDATA ////
    }
}



// Use a single GBT tree to produce a prediction from a single data sample
// of an arbitrary data set. This wouldn't have to be the train set. // IF we use another treeData_ !!
num_t GBT::predictSampleByTree(size_t sampleIdx, size_t treeIdx)
{
  // root of current tree
  Node *currentNode = rootNodes_[treeIdx]; // OLD: &forest_[treeIdx][0];
  num_t value;
  // traverse to the leaf of the tree
  while(currentNode->hasChildren())
    {
      size_t featureIdx = currentNode->getSplitter();
      treeData_->getFeatureData(featureIdx, sampleIdx, value);
      currentNode = currentNode->percolateData(value);
    }
  // now currentNode points to a leaf node: get the prediction
  return currentNode->getPrediction();
}



// Use a single GBT tree to produce predictions for a whole training data set
void GBT::predictDatasetByTree(size_t treeIdx, vector<num_t>& curPrediction)
{
  size_t nSamples = treeData_->nSamples();
  // predict for all samples
  for(size_t i = 0; i < nSamples; ++i)
    {
      curPrediction[i] = GBT::predictSampleByTree(i, treeIdx);
      // cout << "Sample " << i << ", prediction " << curPrediction[i]  << endl;
    }
}



void GBT::growTree(size_t treeIdx)
{
  nLeavesCounter_ = 0; // Counts the number of leaves in the growing tree
  //size_t nSamples = treeData_->nSamples();

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // // ADDED A BOOTSTRAP FUNCTION HERE -- MAYBE IT'S OF USE
  // // NOTE: bootstrapFromRealSamples() generates sample with no NaN's!
  //
  vector<size_t> bootstrapIcs;
  vector<size_t> oobIcs;
  bool withReplacement = false;
  //num_t sampleSize = 0.5; 
  //
  // // 50% of the real samples the target has will be resampled without replacement
  // // Samples left out from bootstrapIcs will be stored in oobIcs
  treeData_->bootstrapFromRealSamples(withReplacement,subSampleSize_,targetIdx_,bootstrapIcs,oobIcs);
  
  if(false)
    {
      vector<num_t> targetData;
      treeData_->getFeatureData(targetIdx_,bootstrapIcs,targetData);
      for(size_t i = 0; i < targetData.size(); ++i)
        {
          assert(!datadefs::isNAN(targetData[i]));
        }

      treeData_->getFeatureData(targetIdx_,oobIcs,targetData);
      for(size_t i = 0; i < targetData.size(); ++i)
        {
          assert(!datadefs::isNAN(targetData[i]));
        }
      cout << "the generated bootstrap sample for tree " << treeIdx << " looks ok" << endl;
    }



  if(false)
    {
      cout << "tree " << treeIdx << "  bootstrap indices [";
      for(size_t i = 0; i < bootstrapIcs.size(); ++i)
        {
          cout << " " << bootstrapIcs[i];
        }
      cout << " ]  oob [";
      for(size_t i = 0; i < oobIcs.size(); ++i)
        {
          cout << " " << oobIcs[i];
        }
      cout << " ]" << endl << endl;
    }


  //////////////////////////////////////////////////////////////////////////////////////////////////

  //Generate a new random  vector for sample indices
  //vector<size_t> sampleIcs(nSamples);
  //for(size_t i = 0; i < nSamples; ++i)
  //  {
  //    sampleIcs[i] = i;
  //  }

  // treeData_->permute(sampleIcs);
  // TODO implement training using only subSampleSize_ samples.

  cout << "Growing tree " << treeIdx << " with a sample:" << endl;

  //size_t rootnode = 0;

  //Start the recursive node splitting from the root node. This will generate the tree.
  // NOTE: COULD sampleIcs BE REPLACED WITH bootstrapIcs?
  GBT::recursiveNodeSplit(treeIdx, rootNodes_[treeIdx], bootstrapIcs);
}



void GBT::SetNodePrediction( size_t treeIdx, Node* node, vector<size_t>& sampleIcs)
{
  // calculate the prediction from leaf node training set
  vector<num_t> data(sampleIcs.size());
  treeData_->getFeatureData(targetIdx_, sampleIcs, data);

  num_t nodePred;
  if (0==numClasses_)
    {
      // numerical target
      size_t nreal;
      datadefs::mean( data, nodePred, nreal); 
    }
  else
    {
      // categorical target, we are predicting one class probability, calculate gamma
      num_t numerator=0.0, denominator=0.0;
      for (size_t i = 0; i < data.size(); ++i)
        {
          num_t abs_data_i = fabs( data[i] );
          denominator += abs_data_i * (1.0 - abs_data_i);
          numerator   += data[i];
        }
      if ( fabs(denominator) <= datadefs::EPS )
        {
          nodePred = datadefs::LOG_OF_MAX_NUM * numerator;
        }
      else
        {
          nodePred = (numClasses_ - 1)*numerator / (numClasses_ * denominator);
        }
    }

  node->setPrediction(nodePred);
  cout << "setting pred=" << nodePred << ", tree="<<treeIdx<<" &node="<< node<<endl;
  
  ++nLeavesCounter_; // finally, indicate that we have added a leaf
}



void GBT::recursiveNodeSplit(size_t treeIdx, Node* node, vector<size_t>& sampleIcs)
{
  size_t nSamples = sampleIcs.size();
  size_t nFeatures = treeData_->nFeatures();

  // This is the GBT termination criterion:
  // Do we exceed the max number of leaves if we split this node? 
  if ( 2+nLeavesCounter_ >= nMaxLeaves_) 
    {
      cout << "tree=" << treeIdx << " &node=" << &node << " nMaxLeaves=" << nMaxLeaves_ << " nLeavesCounter=" << nLeavesCounter_ << endl;
      // Leaf node
      GBT::SetNodePrediction(treeIdx, node, sampleIcs);
      return;
    }
 
  // not enough samples to split?
  if (sampleIcs.size() <= nodeSize_)
    {
      cout << "tree=" << treeIdx << "  &node=" << &node << ". Only " << sampleIcs.size() << " samples!" << endl;
      // Leaf node
      GBT::SetNodePrediction(treeIdx, node, sampleIcs);
      return;
    }
  cout << "tree=" << treeIdx << "  &node=" << &node << " finding split:"<< endl;



  // Find the best target split
  vector<num_t> targetData;
  const bool isTargetNumerical = treeData_->isFeatureNumerical(targetIdx_);
  treeData_->getFeatureData(targetIdx_,sampleIcs,targetData);

  vector<size_t> sampleIcs_left, sampleIcs_right;
  num_t splitValue;
  set<num_t> values_left;


  // Target is ALWAYS numerical with GBT!
  assert( isTargetNumerical == true );
  treeData_->numericalFeatureSplit(targetData,isTargetNumerical,targetData,nodeSize_,sampleIcs_left,sampleIcs_right,splitValue);
  cout << "Target split with itself." << endl;


  // Find the best feature split corresponding to the target split
  num_t bestFitness = 0.0;
  size_t bestFeatureIdx = nFeatures;
  vector<num_t> featureData;


  // TODO Instead of going through all features, sample features here
  // using sampling weights proportional to current importances
  for (size_t featureIdx = 0; featureIdx < nFeatures; ++featureIdx)
    {
      if (featureIdx == targetIdx_)
        {
          continue;
        }
      bool isFeatureNumerical = treeData_->isFeatureNumerical(featureIdx);
      treeData_->getFeatureData(featureIdx,sampleIcs,featureData);
      num_t fitness = treeData_->splitFitness(featureData,isFeatureNumerical,nodeSize_,sampleIcs_left,sampleIcs_right);

      if (fitness > bestFitness)
        {
          bestFitness = fitness;
          bestFeatureIdx = featureIdx;
          if (fabs(bestFitness - 1.0) < datadefs::EPS)
            {
              //cout << "Maximum fitness reached for splitter " << bestfeatureidx << ", stopped searching" << endl;
              break;
            }
        }
    }

  if (bestFeatureIdx == nFeatures)
    {
      cout << "No splitter found, quitting." << endl;
      // Leaf node
      GBT::SetNodePrediction(treeIdx, node, sampleIcs);
      return;
    }
  cout << "Best splitter feature is " << bestFeatureIdx << " with fitness of " << bestFitness << endl;


  treeData_->getFeatureData(bestFeatureIdx,sampleIcs,featureData);

  size_t nRealSamples = 0;
  for(size_t i = 0; i < nSamples; ++i)
    {
      if(!datadefs::isNAN(featureData[i]))
        {
          featureData[nRealSamples] = featureData[i];
          targetData[nRealSamples] = targetData[i];
          ++nRealSamples;
        }
    }
  featureData.resize(nRealSamples);
  targetData.resize(nRealSamples);

  if (nRealSamples < 2*nodeSize_)
    {
      cout << "Splitter has too few non-missing values, quitting" << endl;
      //TODO This needs to be fixed such that one of the surrogates will determine the split instead
      GBT::SetNodePrediction(treeIdx, node, sampleIcs);
      return;
    }


  bool isBestFeatureNumerical = treeData_->isFeatureNumerical(bestFeatureIdx);
  if (isBestFeatureNumerical)
    {
      treeData_->numericalFeatureSplit(targetData,isTargetNumerical,featureData,nodeSize_,sampleIcs_left,sampleIcs_right,splitValue);
    }
  else
    {
      treeData_->categoricalFeatureSplit(targetData,isTargetNumerical,featureData,sampleIcs_left,sampleIcs_right,values_left);
    }

  assert(sampleIcs_left.size() + sampleIcs_right.size() == nRealSamples);

  if (sampleIcs_left.size() < nodeSize_ || sampleIcs_right.size() < nodeSize_)
    {
      // Leaf node
      GBT::SetNodePrediction(treeIdx, node, sampleIcs);
      return;
    }


  if (isBestFeatureNumerical)
    {
      node->setSplitter(bestFeatureIdx,splitValue);
    }
  else
    {
      node->setSplitter(bestFeatureIdx,values_left);
    }


  //WILL BECOME DEPRECATED
  for (size_t i = 0; i < sampleIcs_left.size(); ++i)
  {
    sampleIcs_left[i] = sampleIcs[sampleIcs_left[i]];
  }
  GBT::recursiveNodeSplit(treeIdx,node->leftChild(),sampleIcs_left);

  //WILL BECOME DEPRECATED
  for (size_t i = 0; i < sampleIcs_right.size(); ++i)
  {
    sampleIcs_right[i] = sampleIcs[sampleIcs_right[i]];
  }
  GBT::recursiveNodeSplit(treeIdx,node->rightChild(),sampleIcs_right);
}
