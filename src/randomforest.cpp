#include "randomforest.hpp"
#include "datadefs.hpp"
#include <cmath>
#include <ctime>
#include <iostream>
//#include <omp.h>

Randomforest::Randomforest(Treedata* treedata, 
			   size_t targetIdx, 
			   size_t nTrees, 
			   size_t mTry, 
			   size_t nodeSize,
			   bool useContrasts,
			   bool isOptimizedNodeSplit):
  targetIdx_(targetIdx),
  treedata_(treedata),
  nTrees_(nTrees),
  rootNodes_(nTrees),
  oobMatrix_(nTrees)
{

  if(useContrasts)
    {
      treedata_->permuteContrasts();
    }

  featuresInForest_.clear();

  //These parameters, and those specified in the Random Forest initiatialization, define the type of the forest generated (an RF) 
  bool sampleWithReplacement = true;
  num_t sampleSizeFraction = 1.0;
  size_t maxNodesToStop = treedata->nSamples();
  size_t minNodeSizeToStop = nodeSize;
  bool isRandomSplit = true;
  size_t nFeaturesForSplit = mTry;
  
  //Allocates memory for the root nodes
  for(size_t treeIdx = 0; treeIdx < nTrees_; ++treeIdx)
    {
      rootNodes_[treeIdx] = new RootNode(sampleWithReplacement,
					 sampleSizeFraction,
					 maxNodesToStop,
					 minNodeSizeToStop,
					 isRandomSplit,
					 nFeaturesForSplit,
					 useContrasts,
					 isOptimizedNodeSplit);
    }

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

/*
  size_t Randomforest::getTarget()
  {
  return(targetIdx_);
  }
*/

void Randomforest::growForest()
{

  size_t nNodes;

  void (*leafPredictionFunction)(const vector<num_t>&, num_t&);

  if(treedata_->isFeatureNumerical(targetIdx_))
    {
      leafPredictionFunction = &datadefs::mean;
    }
  else
    {
      leafPredictionFunction = &datadefs::mode;
    }

  //#pragma omp parallel for
  for(size_t treeIdx = 0; treeIdx < nTrees_; ++treeIdx)
    {
      rootNodes_[treeIdx]->growTree(treedata_,
				    targetIdx_,
				    leafPredictionFunction,
				    oobMatrix_[treeIdx],
				    featuresInForest_[treeIdx],
				    nNodes);
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
      num_t leafPrediction = it->first->getLeafTrainPrediction();
      num_t leafImpurity = 0;
      size_t nSamplesInLeaf = targetData.size();
      

      if(isTargetNumerical)
	{
	  for(size_t i = 0; i < nSamplesInLeaf; ++i)
	    {
	      leafImpurity += pow(leafPrediction - targetData[i],2);
	    }
	  //leafImpurity /= nSamplesInLeaf;
	}
      else
	{
	  for(size_t i = 0; i < nSamplesInLeaf; ++i)
	    {
	      if(leafPrediction != targetData[i])
		{
		  ++leafImpurity;
		}
	    }
	  //leafImpurity *= leafImpurity;
	}

      //num_t impurity_leaf;
      //size_t nRealSamples;
      //treedata_->impurity(targetData,isTargetNumerical,impurity_leaf,nRealSamples);
      //size_t n(it->second.size());
      n_tot += nSamplesInLeaf;
      impurity += nSamplesInLeaf * leafImpurity;
    }

  if(n_tot > 0)
    {
      impurity /= n_tot;
    }
}
