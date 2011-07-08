#include "randomforest.hpp"
#include "datadefs.hpp"
#include <cmath>
#include <ctime>
#include <iostream>
//#include <omp.h>

Randomforest::Randomforest(Treedata* treedata, size_t targetIdx, size_t nTrees, size_t mTry, size_t nodeSize):
  targetIdx_(targetIdx),
  treedata_(treedata),
  nTrees_(nTrees),
  mTry_(mTry),
  nodeSize_(nodeSize),
  rootNodes_(nTrees),
  oobMatrix_(nTrees)
{

  featuresInForest_.clear();

  //Allocates memory for the root nodes
  for(size_t treeIdx = 0; treeIdx < nTrees_; ++treeIdx)
    {
      rootNodes_[treeIdx] = new Node;
      //oobMatrix_[treeIdx].clear();
    }

  //cout << "foo1" << endl;

  //Before analysis, we'll permute the contrast features
  treedata_->permuteContrasts();

  //cout << "foo2" << endl;

  //Let's grow the forest
  Randomforest::growForest();

}

Randomforest::~Randomforest()
{

  //cout << "Destructor invoked" << endl;

  //Delete the rootnodes, which will initiate a cascade destroying all the child nodes, thus, the whole trees
  for(size_t treeIdx = 0; treeIdx < nTrees_; ++treeIdx)
    {
      //cout << "deleting rootnode of tree " << treeIdx << endl;
      //cout << " " << forest_[treeIdx].size();
      delete rootNodes_[treeIdx];
      //forest_[treeIdx][0] = NULL;
    }

  //cout << "done" << endl;
} 

size_t Randomforest::getTarget()
{
  return(targetIdx_);
}

void Randomforest::growForest()
{

  //#pragma omp parallel for
  for(size_t treeIdx = 0; treeIdx < nTrees_; ++treeIdx)
    {
      Randomforest::growTree(treeIdx);
    }
    
}

void Randomforest::growTree(size_t treeIdx)
{
  //Generate the vector for bootstrap indices
  vector<size_t> bootstrapIcs;
  bool withReplacement = true;
  num_t sampleSize = 1.0;

  //Generate bootstrap indices and oob-indices
  treedata_->bootstrapFromRealSamples(withReplacement, sampleSize, targetIdx_, bootstrapIcs, oobMatrix_[treeIdx]);

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

  
  //Start the recursive node splitting from the root node. This will generate the tree.
  Randomforest::recursiveNodeSplit(treeIdx,rootNodes_[treeIdx],bootstrapIcs);

}

void Randomforest::recursiveNodeSplit(const size_t treeIdx, Node* node, const vector<size_t>& sampleIcs)
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

  vector<num_t> fitness(mTry_, 0.0);
  num_t bestFitness = 0.0;
  size_t bestFeatureIdx = mTry_;

  //const size_t halfWay = mTry_ / 2;

  //vector<num_t> featureData;
  //omp_lock_t lock;
  //omp_init_lock(&lock);
  //#pragma omp parallel for
  for(int i = 0; i < static_cast<int>(mTry_); ++i)
    {
      //printf("%i\n",i);
      vector<num_t> featureData;
      size_t featureIdx = mTrySample[i];
      bool isFeatureNumerical = treedata_->isFeatureNumerical(featureIdx);

      //Neither the real nor the contrast feature can appear in the tree as splitter
      if(featureIdx == targetIdx_)
        {
          continue;
        }

      //omp_set_lock(&lock);
      treedata_->getFeatureData(featureIdx,sampleIcs,featureData);
      //omp_unset_lock(&lock);

      //omp_set_lock(&lock);
      fitness[i] = treedata_->splitFitness(featureData,isFeatureNumerical,nodeSize_,sampleIcs_left,sampleIcs_right);
      //omp_unset_lock(&lock);

      //omp_set_lock(&lock);
      
      //if(fitness > bestFitness)
      //	{
      //	  bestFitness = fitness;
      //	  bestFeatureIdx = featureIdx;
	  //if(fabs(bestFitness - 1.0) < datadefs::EPS)
	  //  {
	  //cout << "Maximum fitness reached for splitter " << bestfeatureidx << ", stopped searching" << endl;
	  //	break;
	  // }
      //	}
      //omp_unset_lock(&lock);
       
    }
  //omp_destroy_lock(&lock);

  for(size_t i = 0; i < mTry_; ++i)
    {
      if(fitness[i] > bestFitness)
	{
	  bestFitness = fitness[i];
	  bestFeatureIdx = mTrySample[i];
	}
    }


  if(bestFeatureIdx == mTry_)
    {
      //cout << "No splitter found, quitting" << endl << endl;
      return;
    }
  
  //cout << "Best splitter feature is " << bestFeatureIdx << " with fitness of " << bestFitness << endl; 

  vector<num_t> featureData;
  treedata_->getFeatureData(bestFeatureIdx,sampleIcs,featureData);
  
  //vector<size_t> NANIcs;
  
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

  if(nRealSamples < 2*nodeSize_)
    {
      //cout << "Splitter has too few non-missing values, quitting" << endl;
      //This needs to be fixed such that one of the surrogates will determine the split instead
      return;
    }

  //size_t nodeIdx_left; 
  //size_t nodeIdx_right;
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
  
  //cout << "setting splitter" << endl;

  if(isBestFeatureNumerical)
    {
      node->setSplitter(bestFeatureIdx,splitValue);
    }
  else
    {
      node->setSplitter(bestFeatureIdx,values_left);
    }

  featuresInForest_[treeIdx].insert(bestFeatureIdx);
  
  //cout << "setting pointers to new nodes" << endl;

  //assert(nNodes_[treeIdx] == forest_[treeIdx].size());

  

  //cout << nodeIdx_left << " " << nodeIdx_right << endl;

  if(sampleIcs_left.size() > 2*nodeSize_)
    {
      //WILL BECOME DEPRECATED
      for(size_t i = 0; i < sampleIcs_left.size(); ++i)
	{
	  sampleIcs_left[i] = sampleIcs[sampleIcs_left[i]];
	}

      Randomforest::recursiveNodeSplit(treeIdx,node->leftChild(),sampleIcs_left);
    }

  if(sampleIcs_right.size() > 2*nodeSize_)
    {
      //WILL BECOME DEPRECATED
      for(size_t i = 0; i < sampleIcs_right.size(); ++i)
        {
          sampleIcs_right[i] = sampleIcs[sampleIcs_right[i]];
        }

      Randomforest::recursiveNodeSplit(treeIdx,node->rightChild(),sampleIcs_right);
    }
  
}


void Randomforest::percolateSampleIcs(Node* rootNode, vector<size_t>& sampleIcs, map<Node*,vector<size_t> >& trainIcs)
{
  
  trainIcs.clear();
  //map<Node*,vector<size_t> > trainics;
  
  for(size_t i = 0; i < sampleIcs.size(); ++i)
    {
      Node* nodep(rootNode);
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

void Randomforest::percolateSampleIcsAtRandom(size_t featureIdx, Node* rootNode, vector<size_t>& sampleIcs, map<Node*,vector<size_t> >& trainIcs)
{

  trainIcs.clear();

  for(size_t i = 0; i < sampleIcs.size(); ++i)
    {
      Node* nodep(rootNode);
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

vector<num_t> Randomforest::featureImportance()
{

  size_t nRealFeatures = treedata_->nFeatures();
  size_t nAllFeatures = 2*nRealFeatures;
  vector<num_t> importance(nAllFeatures);
  size_t nOobSamples = 0;
  size_t nContrastsInForest = 0;

  for(size_t i = 0; i < nAllFeatures; ++i)
    {
      importance[i] = 0;
    }

  for(map<size_t, set<size_t> >::const_iterator tit(featuresInForest_.begin()); tit != featuresInForest_.end(); ++tit)
    {
      size_t treeIdx = tit->first;
      
      size_t nNewOobSamples = oobMatrix_[treeIdx].size();
      nOobSamples += nNewOobSamples;
      //Node* rootNode(forest_[treeIdx][0]);
      map<Node*,vector<size_t> > trainIcs;
      Randomforest::percolateSampleIcs(rootNodes_[treeIdx],oobMatrix_[treeIdx],trainIcs);
      num_t treeImpurity;
      Randomforest::treeImpurity(trainIcs,treeImpurity);
      //cout << "#nodes_with_train_samples=" << trainics.size() << endl;

      for(set<size_t>::const_iterator fit(tit->second.begin()); fit != tit->second.end(); ++fit)
	{
	  size_t featureIdx = *fit;
	  
	  if(featureIdx >= nRealFeatures)
	    {
	      ++nContrastsInForest;
	    }

	  Randomforest::percolateSampleIcsAtRandom(featureIdx,rootNodes_[treeIdx],oobMatrix_[treeIdx],trainIcs);
	  num_t permutedTreeImpurity;
	  Randomforest::treeImpurity(trainIcs,permutedTreeImpurity);
	  if(fabs(treeImpurity) > datadefs::EPS)
	    {
	      importance[featureIdx] += nNewOobSamples * (permutedTreeImpurity - treeImpurity) / treeImpurity;
	    }
	    
	}
      
    }

  /*
    for(size_t treeIdx = 0; treeIdx < nTrees_; ++treeIdx)
    {
    size_t nNewOobSamples = oobMatrix_[treeIdx].size();
    nOobSamples += nNewOobSamples;
    //Node* rootNode(forest_[treeIdx][0]);
    map<Node*,vector<size_t> > trainIcs;
    Randomforest::percolateSampleIcs(rootNodes_[treeIdx],oobMatrix_[treeIdx],trainIcs);
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
    Randomforest::percolateSampleIcsAtRandom(featureIdx,rootNodes_[treeIdx],oobMatrix_[treeIdx],trainIcs);
    num_t permutedTreeImpurity;
    Randomforest::treeImpurity(trainIcs,permutedTreeImpurity);
    if(fabs(treeImpurity) > datadefs::EPS)
    {
    importance[featureIdx] += nNewOobSamples * (permutedTreeImpurity - treeImpurity) / treeImpurity;
    } 
    }
    }
    }
  */

  size_t nNodesInForest = Randomforest::nNodes();
  
  for(size_t featureIdx = 0; featureIdx < nAllFeatures; ++featureIdx)
    {
      importance[featureIdx] *= 100.0*nTrees_/nNodesInForest;//1.0*nNodesInForest/nTrees_; //nContrastsInForest
      importance[featureIdx] /= nOobSamples; //nRealFeatures
    }

  return(importance);

}

size_t Randomforest::nNodes()
{
  size_t nNodes = 0;
  for(size_t treeIdx = 0; treeIdx < nTrees_; ++treeIdx)
    {
      nNodes += rootNodes_[treeIdx]->nNodes();
    }
  
  return(nNodes);
}

vector<size_t> Randomforest::featureFrequency()
{
  assert(false);
  vector<size_t> frequency(0);
  return(frequency);
}




/*
  bool Randomforest::isFeatureInTree(size_t featureIdx, size_t treeIdx)
  {
  set<size_t>::const_iterator it(featuresInForest_[featureIdx].find(featureIdx));
  if(it == featuresInForest_[treeIdx].end())
  {
  return(false);
  }
  else
  {
  return(true);
  }
  }
*/





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
