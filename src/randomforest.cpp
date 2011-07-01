#include "randomforest.hpp"
#include "datadefs.hpp"
#include <cmath>
#include <ctime>
#include <iostream>

Randomforest::Randomforest(Treedata* treedata, size_t targetIdx, size_t nTrees, size_t mTry, size_t nodeSize):
  targetIdx_(targetIdx),
  treedata_(treedata),
  nTrees_(nTrees),
  mTry_(mTry),
  nodeSize_(nodeSize),
  oobMatrix_(nTrees)
{

  //CHECK THE PURPOSE OF THIS FUNCTION
  Randomforest::setTarget(targetIdx_);

  //cout << forest_.size() << " trees and " << forest_[0].size() << " max nodes per tree initialized." << endl;

}

Randomforest::~Randomforest()
{

} 

void Randomforest::initialize()
{

  size_t nRealSamples = treedata_->nRealSamples(targetIdx_);
  vector<num_t> targetData;
  treedata_->getFeatureData(targetIdx_,targetData);
  datadefs::countRealValues(targetData,nRealSamples);

  //First we count the theoretical maximum number of nodes per tree.
  //Because each leaf must contain at least nodesize amount of data points, nmaxleaves is
  int nMaxLeaves = int(ceil(float(nRealSamples)/nodeSize_));
  //The upper bound for depth of the tree is log2(nmaxleaves)=log10(nmaxleaves)/log10(2.0):
  int maxDepth = int(ceil(log10(float(nMaxLeaves))/log10(2.0)));
  //Thus, the number of nodes in a complete binary tree of depth maxdepth, there are
  int nMaxNodes = int(pow(2.0,maxDepth+1)); //In reality it's nmaxnodes-1 but this way we'll get a power of two which is supposed to be faster :)

  nNodes_.clear();
  forest_.clear();

  nNodes_.resize(nTrees_);
  forest_.resize(nTrees_);
  for(size_t treeIdx = 0; treeIdx < nTrees_; ++treeIdx)
    {
      forest_[treeIdx].resize(nMaxNodes);
      nNodes_[treeIdx] = 1;
      oobMatrix_[treeIdx].clear();
    }

  treedata_->permuteContrasts();
}
  
void Randomforest::setTarget(size_t targetIdx)
{
  targetIdx_ = targetIdx;
  Randomforest::initialize();
}

size_t Randomforest::getTarget()
{
  return(targetIdx_);
}

void Randomforest::learn(const size_t nPerms, vector<num_t>& pValues, vector<num_t>& importanceValues)
{
  assert(nPerms > 5);
  vector<vector<num_t> > importanceMat(nPerms);
  
  clock_t time_start(clock());
  size_t nNodesInAllForests = 0;
  for(size_t permIdx = 0; permIdx < nPerms; ++permIdx)
    {
      cout << "  RF " << permIdx + 1 << ": ";
      Randomforest::initialize();
      treedata_->permuteContrasts();
      size_t nNodesInForest = 0;
      
      for(size_t treeIdx = 0; treeIdx < nTrees_; ++treeIdx)
	{
	  Randomforest::growTree(treeIdx);
	  nNodesInForest += nNodes_[treeIdx];
	}
      
      nNodesInAllForests += nNodesInForest;
      
      Randomforest::calculateFeatureImportance(importanceMat[permIdx]); //REWORK!
      
      cout << nNodesInForest << " nodes (avg. " << 1.0*nNodesInForest / nTrees_ << " nodes / tree)" << endl;
    }

  num_t time_diff = 1.0*(clock() - time_start) / CLOCKS_PER_SEC;
  cout << nPerms << " RFs, " << nPerms*nTrees_ << " trees, and " << nNodesInAllForests 
       << " nodes generated in " << time_diff << " seconds (" << 1.0*nNodesInAllForests / time_diff 
       << " nodes per second)" << endl;

  size_t nFeatures = treedata_->nFeatures();
  pValues.resize(nFeatures);
  
  for(size_t featureIdx = 0; featureIdx < nFeatures; ++featureIdx)
    {

      size_t nRealSamples;
      vector<num_t> fSample(nPerms);
      vector<num_t> cSample(nPerms);
      for(size_t permIdx = 0; permIdx < nPerms; ++permIdx)
	{
	  fSample[permIdx] = importanceMat[permIdx][featureIdx];
	  cSample[permIdx] = importanceMat[permIdx][featureIdx + nFeatures];
	}
      datadefs::utest(fSample,cSample,pValues[featureIdx]);
      datadefs::mean(fSample,importanceValues[featureIdx],nRealSamples);
    }
  

  //cout << "done" << endl;
  
}

void Randomforest::growTree(size_t treeIdx)
{
  //Generate the vector for bootstrap indices
  vector<size_t> bootstrapIcs;

  //Generate bootstrap indices and oob-indices
  treedata_->bootstrapFromRealSamples(targetIdx_,bootstrapIcs,oobMatrix_[treeIdx]);

  //This is to check that the bootstrap sample doesn't contain any missing values (it shouldn't!)
  if(false)
    {
      vector<num_t> targetData;
      treedata_->getFeatureData(targetIdx_,bootstrapIcs,targetData);
      for(size_t i = 0; i < targetData.size(); ++i)
	{
	  assert(!datadefs::isNAN(targetData[i]));
	}
      
      treedata_->getFeatureData(targetIdx_,oobMatrix_[treeIdx],targetData);
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
      for(size_t i = 0; i < oobMatrix_[treeIdx].size(); ++i)
	{
	  cout << " " << oobMatrix_[treeIdx][i];
	}
      cout << " ]" << endl << endl;
    }

  size_t rootNode = 0;
  
  //Start the recursive node splitting from the root node. This will generate the tree.
  Randomforest::recursiveNodeSplit(treeIdx,rootNode,bootstrapIcs);

}

void Randomforest::recursiveNodeSplit(const size_t treeIdx, const size_t nodeIdx, const vector<size_t>& sampleIcs)
{

  size_t nSamples = sampleIcs.size();

  if(nSamples < 2*nodeSize_)
    {
      //cout << "Too few samples to start with, quitting" << endl;
      return;
    }

  //size_t nfeatures = treedata_->nfeatures();
  //Create mtry randomly selected feature indices to determine the split
  vector<size_t> mTrySample(mTry_);
  
  const size_t nFeatures = treedata_->nFeatures();

  for(size_t i = 0; i < mTry_; ++i)
    {
      treedata_->getRandomIndex(2*nFeatures,mTrySample[i]);
    }

  //const size_t targetIdx = treedata_->getTarget();
  vector<num_t> targetData;
  const bool isTargetNumerical = treedata_->isFeatureNumerical(targetIdx_);
  treedata_->getFeatureData(targetIdx_,sampleIcs,targetData);

  //assert(!datadefs::isNAN(targetData));

  vector<size_t> sampleIcs_left,sampleIcs_right;
  num_t splitValue;
  set<num_t> values_left;

  if(isTargetNumerical)
    {
      treedata_->numericalFeatureSplit(targetData,isTargetNumerical,targetData,nodeSize_,sampleIcs_left,sampleIcs_right,splitValue);
    }
  else
    {
      treedata_->categoricalFeatureSplit(targetData,isTargetNumerical,targetData,sampleIcs_left,sampleIcs_right,values_left);
    }

  //cout << "Target splitted with itself." << endl;

  num_t bestFitness = 0.0;
  size_t bestFeatureIdx = mTry_;

  //const size_t halfWay = mTry_ / 2;

  vector<num_t> featureData;
  for(size_t i = 0; i < mTry_; ++i)
    {
      size_t featureIdx = mTrySample[i];
      bool isFeatureNumerical = treedata_->isFeatureNumerical(featureIdx);

      //Neither the real nor the contrast feature can appear in the tree as splitter
      if(featureIdx == targetIdx_)
        {
          continue;
        }

      //First half of mtry are real features, the other half contrast features
      //if(i < halfWay)
      //	{
      treedata_->getFeatureData(featureIdx,sampleIcs,featureData);
      //	}
      //else
      //	{
      //treedata_->getContrastData(featureIdx,sampleIcs,featureData);
      //	}

      num_t fitness = treedata_->splitFitness(featureData,isFeatureNumerical,nodeSize_,sampleIcs_left,sampleIcs_right);

      if(fitness > bestFitness)
	{
	  bestFitness = fitness;
	  bestFeatureIdx = featureIdx;
	  if(fabs(bestFitness - 1.0) < datadefs::EPS)
	    {
	      //cout << "Maximum fitness reached for splitter " << bestfeatureidx << ", stopped searching" << endl;
	      break;
	    }
	}
    }
  
  if(bestFeatureIdx == mTry_)
    {
      //cout << "No splitter found, quitting" << endl << endl;
      return;
    }
  
  //cout << "Best splitter feature is " << bestFeatureIdx << " with fitness of " << bestFitness << endl; 

  treedata_->getFeatureData(bestFeatureIdx,sampleIcs,featureData);
  
  vector<size_t> NANIcs;
  
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

  //cout << "Splitter " << bestFeatureIdx << " has " << NANIcs.size() << " missing values, which will be omitted in splitting" << endl;

  /*
    for(int i = NANIcs.size()-1; i >= 0; --i)
    {
    size_t removeIdx = NANIcs[i];
    //cout << "removing sample " << removeIdx << ": " << targetData[removeIdx] << " and " << featureData[removeIdx] << endl;
    targetData.erase(targetData.begin() + removeIdx);
    featureData.erase(featureData.begin() + removeIdx);
    //sampleIcs.erase(sampleIcs.begin() + i);
    }
  */
  //n_tot = targetData.size();

  //assert(!datadefs::isNAN(featureData));
  //assert(!datadefs::isNAN(targetData));

  if(nRealSamples < 2*nodeSize_)
    {
      //cout << "Splitter has too few non-missing values, quitting" << endl;
      //This needs to be fixed such that one of the surrogates will determine the split instead
      return;
    }

  size_t nodeIdx_left; 
  size_t nodeIdx_right;
  bool isBestFeatureNumerical = treedata_->isFeatureNumerical(bestFeatureIdx);

  if(isBestFeatureNumerical)
    {
      treedata_->numericalFeatureSplit(targetData,isTargetNumerical,featureData,nodeSize_,sampleIcs_left,sampleIcs_right,splitValue);
    }
  else
    {
      treedata_->categoricalFeatureSplit(targetData,isTargetNumerical,featureData,sampleIcs_left,sampleIcs_right,values_left);
    }

  //assert(sampleIcs.size() == n_tot);
  assert(sampleIcs_left.size() + sampleIcs_right.size() == nRealSamples);

  if(sampleIcs_left.size() < nodeSize_ || sampleIcs_right.size() < nodeSize_)
    {
      //cout << "Split was unsuccessful, quitting" << endl;
      return;
    }
  
  nodeIdx_left = nNodes_[treeIdx]++;
  nodeIdx_right = nNodes_[treeIdx]++;
  
  if(isBestFeatureNumerical)
    {
      forest_[treeIdx][nodeIdx].setSplitter(bestFeatureIdx,splitValue,forest_[treeIdx][nodeIdx_left],forest_[treeIdx][nodeIdx_right]);
    }
  else
    {
      forest_[treeIdx][nodeIdx].setSplitter(bestFeatureIdx,values_left,forest_[treeIdx][nodeIdx_left],forest_[treeIdx][nodeIdx_right]);
    }
  
  if(sampleIcs_left.size() > 2*nodeSize_)
    {
      //WILL BECOME DEPRECATED
      for(size_t i = 0; i < sampleIcs_left.size(); ++i)
	{
	  sampleIcs_left[i] = sampleIcs[sampleIcs_left[i]];
	}

      Randomforest::recursiveNodeSplit(treeIdx,nodeIdx_left,sampleIcs_left);
    }

  if(sampleIcs_right.size() > 2*nodeSize_)
    {
      //WILL BECOME DEPRECATED
      for(size_t i = 0; i < sampleIcs_right.size(); ++i)
        {
          sampleIcs_right[i] = sampleIcs[sampleIcs_right[i]];
        }

      Randomforest::recursiveNodeSplit(treeIdx,nodeIdx_right,sampleIcs_right);
    }
  
}


void Randomforest::percolateSampleIcs(Node& rootNode, vector<size_t>& sampleIcs, map<Node*,vector<size_t> >& trainIcs)
{
  
  trainIcs.clear();
  //map<Node*,vector<size_t> > trainics;
  
  for(size_t i = 0; i < sampleIcs.size(); ++i)
    {
      Node* nodep(&rootNode);
      size_t sampleIdx = sampleIcs[i];
      Randomforest::percolateSampleIdx(sampleIdx,&nodep);
      map<Node*,vector<size_t> >::iterator it(trainIcs.find(nodep));
      if(it == trainIcs.end())
	{
	  Node* foop(nodep);
	  vector<size_t> foo(1);
	  foo[0] = sampleIdx;
	  trainIcs.insert(pair<Node*,vector<size_t> >(foop,foo));
	}
      else
	{
	  trainIcs[it->first].push_back(sampleIdx);
	}
      
    }
  
  
  if(false)
    {
      cout << "Train samples percolated accordingly:" << endl;
      size_t iter = 0;
      for(map<Node*,vector<size_t> >::const_iterator it(trainIcs.begin()); it != trainIcs.end(); ++it, ++iter)
	{
	  cout << "leaf node " << iter << ":"; 
	  for(size_t i = 0; i < it->second.size(); ++i)
	    {
	      cout << " " << it->second[i];
	    }
	  cout << endl;
	}
    }
}

void Randomforest::percolateSampleIcsAtRandom(size_t featureIdx, Node& rootNode, vector<size_t>& sampleIcs, map<Node*,vector<size_t> >& trainIcs)
{

  trainIcs.clear();

  for(size_t i = 0; i < sampleIcs.size(); ++i)
    {
      Node* nodep(&rootNode);
      size_t sampleIdx = sampleIcs[i];
      Randomforest::percolateSampleIdxAtRandom(featureIdx,sampleIdx,&nodep);
      map<Node*,vector<size_t> >::iterator it(trainIcs.find(nodep));
      if(it == trainIcs.end())
        {
          Node* foop(nodep);
          vector<size_t> foo(1);
          foo[0] = sampleIdx;
          trainIcs.insert(pair<Node*,vector<size_t> >(foop,foo));
        }
      else
        {
          trainIcs[it->first].push_back(sampleIdx);
        }

    }
}

void Randomforest::percolateSampleIdx(size_t sampleIdx, Node** nodep)
{
  while((*nodep)->hasChildren())
    {
      int featureIdxNew((*nodep)->getSplitter());
      num_t value;
      treedata_->getFeatureData(featureIdxNew,sampleIdx,value);
      *nodep = (*nodep)->percolateData(value);
    }
}

void Randomforest::percolateSampleIdxAtRandom(size_t featureIdx, size_t sampleIdx, Node** nodep)
{
  while((*nodep)->hasChildren())
    {
      size_t featureIdxNew = (*nodep)->getSplitter();
      num_t value = datadefs::NUM_NAN;
      if(featureIdx == featureIdxNew)
	{
	  while(datadefs::isNAN(value))
	    {
	      treedata_->getRandomData(featureIdxNew,value);
	    }
	}
      else
	{
	  treedata_->getFeatureData(featureIdxNew,sampleIdx,value);
	}
      *nodep = (*nodep)->percolateData(value);
    }
}

//NEEDS REWORKING
void Randomforest::calculateFeatureImportance(vector<num_t>& importance)
{

  size_t nRealFeatures = treedata_->nFeatures() ;
  size_t nAllFeatures = 2*nRealFeatures ;
  importance.resize(nAllFeatures);
  size_t nOobSamples = 0;

  size_t nContrastsInForest = 0;

  for(size_t i = 0; i < nAllFeatures; ++i)
    {
      importance[i] = 0;
    }

  for(size_t treeIdx = 0; treeIdx < nTrees_; ++treeIdx)
    {
      size_t nNewOobSamples = oobMatrix_[treeIdx].size();
      nOobSamples += nNewOobSamples;
      Node rootNode(forest_[treeIdx][0]);
      map<Node*,vector<size_t> > trainIcs;
      Randomforest::percolateSampleIcs(rootNode,oobMatrix_[treeIdx],trainIcs);
      num_t treeImpurity;
      Randomforest::treeImpurity(trainIcs,treeImpurity);
      //cout << "#nodes_with_train_samples=" << trainics.size() << endl;  
      for(size_t featureIdx = 0; featureIdx < nAllFeatures; ++featureIdx)
	{
	  if(Randomforest::isFeatureInTree(featureIdx,treeIdx))
	    {
	      if(featureIdx >= nRealFeatures)
		{
		  ++nContrastsInForest;
		}
	      Randomforest::percolateSampleIcsAtRandom(featureIdx,rootNode,oobMatrix_[treeIdx],trainIcs);
	      num_t permutedTreeImpurity;
	      Randomforest::treeImpurity(trainIcs,permutedTreeImpurity);
	      if(fabs(treeImpurity) > datadefs::EPS)
		{
		  importance[featureIdx] += nNewOobSamples * (permutedTreeImpurity - treeImpurity) / treeImpurity;
		} 
	    }
	}
    }

  size_t nNodesInForest = 0;
  for(size_t i = 0; i < nNodes_.size(); ++i)
    {
      nNodesInForest += nNodes_[i];
    }

  for(size_t featureIdx = 0; featureIdx < nAllFeatures; ++featureIdx)
    {
      importance[featureIdx] *= 1.0*nNodesInForest/nTrees_; //nContrastsInForest
      importance[featureIdx] /= nOobSamples; //nRealFeatures
    }

}

bool Randomforest::isFeatureInTree(size_t featureIdx, size_t treeIdx)
{
  for(size_t nodeIdx = 0; nodeIdx < nNodes_[treeIdx]; ++nodeIdx)
    {
      if(forest_[treeIdx][nodeIdx].hasChildren())
      	{
      if(featureIdx == forest_[treeIdx][nodeIdx].getSplitter())
	{
	  return(true);
	}
	}
    }
  return(false);
}

void Randomforest::treeImpurity(map<Node*,vector<size_t> >& trainIcs, num_t& impurity)
{

  impurity = 0.0;
  size_t n_tot = 0;
  
  //size_t targetIdx = treedata_->getTarget();
  bool isTargetNumerical = treedata_->isFeatureNumerical(targetIdx_);

  for(map<Node*,vector<size_t> >::iterator it(trainIcs.begin()); it != trainIcs.end(); ++it)
    {

      vector<num_t> targetData;
      treedata_->getFeatureData(targetIdx_,it->second,targetData);

      num_t impurity_leaf;
      size_t nRealSamples;
      treedata_->impurity(targetData,isTargetNumerical,impurity_leaf,nRealSamples);
      //size_t n(it->second.size());
      n_tot += nRealSamples;
      impurity += nRealSamples * impurity_leaf;
    }

  if(n_tot > 0)
    {
      impurity /= n_tot;
    }
}
