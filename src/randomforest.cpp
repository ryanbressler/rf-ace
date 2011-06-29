#include "randomforest.hpp"
#include "datadefs.hpp"
#include <cmath>
#include <ctime>
#include <iostream>

Randomforest::Randomforest(Treedata* treedata, size_t ntrees, size_t mtry, size_t nodesize):
  treedata_(treedata),
  ntrees_(ntrees),
  mtry_(mtry),
  nodesize_(nodesize),
  oobmatrix_(ntrees)
{

  size_t defaulttargetidx = 0;
  Randomforest::selectTarget(defaulttargetidx);

  //cout << forest_.size() << " trees and " << forest_[0].size() << " max nodes per tree initialized." << endl;

}

Randomforest::~Randomforest()
{

} 

void Randomforest::initializeForest()
{

  //size_t nsamples = treedata_->nsamples();
  size_t nsamples = treedata_->nrealsamples();

  //First we count the theoretical maximum number of nodes per tree.
  //Because each leaf must contain at least nodesize amount of data points, nmaxleaves is
  int nmaxleaves = int(ceil(float(nsamples)/nodesize_));
  //The upper bound for depth of the tree is log2(nmaxleaves)=log10(nmaxleaves)/log10(2.0):
  int maxdepth = int(ceil(log10(float(nmaxleaves))/log10(2.0)));
  //Thus, the number of nodes in a complete binary tree of depth maxdepth, there are
  int nmaxnodes = int(pow(2.0,maxdepth+1)); //In reality it's nmaxnodes-1 but this way we'll get a power of two which is supposed to be faster :)

  nnodes_.clear();
  forest_.clear();

  nnodes_.resize(ntrees_);
  forest_.resize(ntrees_);
  for(size_t i = 0; i < ntrees_; ++i)
    {
      forest_[i].resize(nmaxnodes);
      nnodes_[i] = 1;
      oobmatrix_[i].clear();
    }
}
  
void Randomforest::selectTarget(size_t targetIdx)
{
  if(treedata_->getTarget() != targetIdx)
    {
      treedata_->selectTarget(targetIdx);
    }

  Randomforest::initializeForest();
}

size_t Randomforest::getTarget()
{
  return(treedata_->getTarget());
}

void Randomforest::growForestEnsemble(const size_t nperms, vector<num_t>& pvalues, vector<num_t>& ivalues)
{
  assert(nperms > 5);
  vector<vector<num_t> > importancemat(nperms);
  //vector<num_t> csample(nperms);

  clock_t time_start(clock());
  size_t nnodesinallforests = 0;
  for(size_t p = 0; p < nperms; ++p)
    {
      cout << "  RF " << p + 1 << ": ";
      Randomforest::initializeForest();
      treedata_->permuteContrasts();
      size_t nnodesinforest = 0;
      for(size_t treeIdx = 0; treeIdx < ntrees_; ++treeIdx)
	{
	  Randomforest::growTree(treeIdx);
	  nnodesinforest += nnodes_[treeIdx];
	}
      nnodesinallforests += nnodesinforest;
      //Randomforest::calculate_importance(alpha,importancemat[p],csample[p]);
      Randomforest::calculate_importance(importancemat[p]);
      cout << nnodesinforest << " nodes (avg. " << 1.0*nnodesinforest / ntrees_ << " nodes / tree)" << endl;
    }

  num_t time_diff = 1.0*(clock() - time_start) / CLOCKS_PER_SEC;
  cout << nperms << " RFs, " << nperms*ntrees_ << " trees, and " << nnodesinallforests 
       << " nodes generated in " << time_diff << " seconds (" << 1.0*nnodesinallforests / time_diff 
       << " nodes per second)" << endl;

  size_t nfeatures = treedata_->nfeatures();
  pvalues.resize(nfeatures);
  
  for(size_t f = 0; f < nfeatures; ++f)
    {

      size_t nreal;
      vector<num_t> fsample(nperms);
      vector<num_t> csample(nperms);
      for(size_t p = 0; p < nperms; ++p)
	{
	  fsample[p] = importancemat[p][f];
	  csample[p] = importancemat[p][f + nfeatures];
	}
      datadefs::utest(fsample,csample,pvalues[f]);
      datadefs::mean(fsample,ivalues[f],nreal);
    }
  

  //cout << "done" << endl;
  
}

void Randomforest::growTree(size_t treeIdx)
{
  //Generate the vector for bootstrap indices
  vector<size_t> bootstrapIcs;
  
  //TODO: Redo generation of bootstrap indices
  //**************************************

  //Generate bootstrap indices and oob-indices
  treedata_->bootstrap(bootstrapIcs,oobmatrix_[treeIdx]);

  //This is to check that the bootstrap sample doesn't contain any missing values (it shouldn't!)
  if(false)
    {
      vector<num_t> targetData;
      treedata_->getFeatureData(treedata_->getTarget(),bootstrapIcs,targetData);
      for(size_t i = 0; i < targetData.size(); ++i)
	{
	  assert(!datadefs::isNAN(targetData[i]));
	}
      
      treedata_->getFeatureData(treedata_->getTarget(),oobmatrix_[treeIdx],targetData);
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
      for(size_t i = 0; i < oobmatrix_[treeIdx].size(); ++i)
	{
	  cout << " " << oobmatrix_[treeIdx][i];
	}
      cout << " ]" << endl << endl;
    }

  size_t rootNode = 0;
  
  //Start the recursive node splitting from the root node. This will generate the tree.
  Randomforest::recursiveNodesplit(treeIdx,rootNode,bootstrapIcs);

}

void Randomforest::recursiveNodesplit(const size_t treeIdx, const size_t nodeIdx, const vector<size_t>& sampleIcs)
{

  size_t n_tot = sampleIcs.size();

  if(n_tot < 2*nodesize_)
    {
      //cout << "Too few samples to start with, quitting" << endl;
      return;
    }

  //size_t nfeatures = treedata_->nfeatures();
  //Create mtry randomly selected feature indices to determine the split
  vector<size_t> mtrysample(mtry_);
  
  const size_t nFeatures = treedata_->nfeatures();

  for(size_t i = 0; i < mtry_; ++i)
    {
      mtrysample[i] = treedata_->sampleRandomIdx(nFeatures);
    }

  const size_t targetIdx = treedata_->getTarget();
  vector<num_t> targetData;
  const bool isTargetNumerical = treedata_->isFeatureNumerical(treedata_->getTarget());
  treedata_->getFeatureData(treedata_->getTarget(),sampleIcs,targetData);

  //assert(!datadefs::isNAN(targetData));

  vector<size_t> sampleIcs_left,sampleIcs_right;
  num_t splitValue;
  set<num_t> values_left;

  if(isTargetNumerical)
    {
      treedata_->numericalFeatureSplit(targetData,isTargetNumerical,targetData,nodesize_,sampleIcs_left,sampleIcs_right,splitValue);
    }
  else
    {
      treedata_->categoricalFeatureSplit(targetData,isTargetNumerical,targetData,sampleIcs_left,sampleIcs_right,values_left);
    }

  //cout << "Target splitted with itself." << endl;

  num_t bestFitness = 0.0;
  size_t bestFeatureIdx = mtry_;

  const size_t halfWay = mtry_ / 2;

  vector<num_t> featureData;
  for(size_t i = 0; i < mtry_; ++i)
    {
      size_t featureIdx = mtrysample[i];
      bool isFeatureNumerical = treedata_->isFeatureNumerical(featureIdx);

      //Neither the real nor the contrast feature can appear in the tree as splitter
      if(featureIdx == targetIdx)
        {
          continue;
        }

      //First half of mtry are real features, the other half contrast features
      if(i < halfWay)
	{
	  treedata_->getFeatureData(featureIdx,sampleIcs,featureData);
	}
      else
	{
	  treedata_->getContrastData(featureIdx,sampleIcs,featureData);
	}

      num_t fitness = treedata_->splitFitness(featureData,isFeatureNumerical,nodesize_,sampleIcs_left,sampleIcs_right);

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
  
  if(bestFeatureIdx == mtry_)
    {
      //cout << "No splitter found, quitting" << endl << endl;
      return;
    }
  
  //cout << "Best splitter feature is " << bestFeatureIdx << " with fitness of " << bestFitness << endl; 

  treedata_->getFeatureData(bestFeatureIdx,sampleIcs,featureData);
  
  vector<size_t> NANIcs;
  datadefs::findNANs(featureData,NANIcs);

  //cout << "Splitter " << bestFeatureIdx << " has " << NANIcs.size() << " missing values, which will be omitted in splitting" << endl;

  for(int i = NANIcs.size()-1; i >= 0; --i)
    {
      size_t removeIdx = NANIcs[i];
      //cout << "removing sample " << removeIdx << ": " << targetData[removeIdx] << " and " << featureData[removeIdx] << endl;
      targetData.erase(targetData.begin() + removeIdx);
      featureData.erase(featureData.begin() + removeIdx);
      //sampleIcs.erase(sampleIcs.begin() + i);
    }
  n_tot = targetData.size();

  //assert(!datadefs::isNAN(featureData));
  //assert(!datadefs::isNAN(targetData));

  if(n_tot < 2*nodesize_)
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
      treedata_->numericalFeatureSplit(targetData,isTargetNumerical,featureData,nodesize_,sampleIcs_left,sampleIcs_right,splitValue);
    }
  else
    {
      treedata_->categoricalFeatureSplit(targetData,isTargetNumerical,featureData,sampleIcs_left,sampleIcs_right,values_left);
    }

  //assert(sampleIcs.size() == n_tot);
  assert(sampleIcs_left.size() + sampleIcs_right.size() == featureData.size());
  assert(sampleIcs_left.size() + sampleIcs_right.size() == targetData.size());

  if(sampleIcs_left.size() < nodesize_ || sampleIcs_right.size() < nodesize_)
    {
      //cout << "Split was unsuccessful, quitting" << endl;
      return;
    }
  
  nodeIdx_left = nnodes_[treeIdx]++;
  nodeIdx_right = nnodes_[treeIdx]++;
  
  if(isBestFeatureNumerical)
    {
      forest_[treeIdx][nodeIdx].set_splitter(bestFeatureIdx,splitValue,forest_[treeIdx][nodeIdx_left],forest_[treeIdx][nodeIdx_right]);
    }
  else
    {
      forest_[treeIdx][nodeIdx].set_splitter(bestFeatureIdx,values_left,forest_[treeIdx][nodeIdx_left],forest_[treeIdx][nodeIdx_right]);
    }
  
  if(sampleIcs_left.size() > 2*nodesize_)
    {
      for(size_t i = 0; i < sampleIcs_left.size(); ++i)
	{
	  sampleIcs_left[i] = sampleIcs[sampleIcs_left[i]];
	}

      Randomforest::recursiveNodesplit(treeIdx,nodeIdx_left,sampleIcs_left);
    }

  if(sampleIcs_right.size() > 2*nodesize_)
    {
      for(size_t i = 0; i < sampleIcs_right.size(); ++i)
        {
          sampleIcs_right[i] = sampleIcs[sampleIcs_right[i]];
        }

      Randomforest::recursiveNodesplit(treeIdx,nodeIdx_right,sampleIcs_right);
    }
  
}


void Randomforest::percolate_sampleics(Node& rootnode, vector<size_t>& sampleics, map<Node*,vector<size_t> >& trainics)
{
  
  trainics.clear();
  //map<Node*,vector<size_t> > trainics;
  
  for(size_t i = 0; i < sampleics.size(); ++i)
    {
      Node* nodep(&rootnode);
      size_t sampleidx(sampleics[i]);
      Randomforest::percolate_sampleidx(sampleidx,&nodep);
      map<Node*,vector<size_t> >::iterator it(trainics.find(nodep));
      if(it == trainics.end())
	{
	  Node* foop(nodep);
	  vector<size_t> foo(1);
	  foo[0] = sampleidx;
	  trainics.insert(pair<Node*,vector<size_t> >(foop,foo));
	}
      else
	{
	  trainics[it->first].push_back(sampleidx);
	}
      
    }
  
  
  if(false)
    {
      cout << "Train samples percolated accordingly:" << endl;
      size_t iter = 0;
      for(map<Node*,vector<size_t> >::const_iterator it(trainics.begin()); it != trainics.end(); ++it, ++iter)
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

void Randomforest::percolate_sampleics_randf(size_t featureidx, Node& rootnode, vector<size_t>& sampleics, map<Node*,vector<size_t> >& trainics)
{

  trainics.clear();

  for(size_t i = 0; i < sampleics.size(); ++i)
    {
      Node* nodep(&rootnode);
      size_t sampleidx(sampleics[i]);
      Randomforest::percolate_sampleidx_randf(featureidx,sampleidx,&nodep);
      map<Node*,vector<size_t> >::iterator it(trainics.find(nodep));
      if(it == trainics.end())
        {
          Node* foop(nodep);
          vector<size_t> foo(1);
          foo[0] = sampleidx;
          trainics.insert(pair<Node*,vector<size_t> >(foop,foo));
        }
      else
        {
          trainics[it->first].push_back(sampleidx);
        }

    }
}

void Randomforest::percolate_sampleidx(size_t sampleidx, Node** nodep)
{
  while((*nodep)->has_children())
    {
      size_t featureidx_new((*nodep)->get_splitter());
      num_t value;
      treedata_->getFeatureData(featureidx_new,sampleidx,value);
      *nodep = (*nodep)->percolate(value);
    }
}

void Randomforest::percolate_sampleidx_randf(size_t featureidx, size_t sampleidx, Node** nodep)
{
  while((*nodep)->has_children())
    {
      size_t featureidx_new((*nodep)->get_splitter());
      num_t value;
      if(featureidx == featureidx_new)
	{
	  value = treedata_->sampleAtRandom(featureidx_new);
	}
      else
	{
	  treedata_->getFeatureData(featureidx_new,sampleidx,value);
	}
      *nodep = (*nodep)->percolate(value);
    }
}


void Randomforest::calculate_importance(vector<num_t>& importance)
{

  size_t nRealFeatures = treedata_->nfeatures() ;
  size_t nAllFeatures = 2*nRealFeatures ;
  importance.resize(nAllFeatures);
  size_t nOobSamples = 0;

  size_t nContrastsInForest = 0;

  for(size_t i = 0; i < nAllFeatures; ++i)
    {
      importance[i] = 0;
    }

  for(size_t treeIdx = 0; treeIdx < ntrees_; ++treeIdx)
    {
      size_t nNewOobSamples = oobmatrix_[treeIdx].size();
      nOobSamples += nNewOobSamples;
      Node rootNode(forest_[treeIdx][0]);
      map<Node*,vector<size_t> > trainIcs;
      Randomforest::percolate_sampleics(rootNode,oobmatrix_[treeIdx],trainIcs);
      num_t treeImpurity;
      Randomforest::tree_impurity(trainIcs,treeImpurity);
      //cout << "#nodes_with_train_samples=" << trainics.size() << endl;  
      for(size_t featureIdx = 0; featureIdx < nAllFeatures; ++featureIdx)
	{
	  if(Randomforest::is_feature_in_tree(featureIdx,treeIdx))
	    {
	      if(featureIdx >= nRealFeatures)
		{
		  ++nContrastsInForest;
		}
	      Randomforest::percolate_sampleics_randf(featureIdx,rootNode,oobmatrix_[treeIdx],trainIcs);
	      num_t permutedTreeImpurity;
	      Randomforest::tree_impurity(trainIcs,permutedTreeImpurity);
	      if(fabs(treeImpurity) > datadefs::EPS)
		{
		  importance[featureIdx] += nNewOobSamples * (permutedTreeImpurity - treeImpurity) / treeImpurity;
		} 
	    }
	}
    }

  size_t nNodesInForest = 0;
  for(size_t i = 0; i < nnodes_.size(); ++i)
    {
      nNodesInForest += nnodes_[i];
    }
  
  for(size_t featureIdx = 0; featureIdx < nAllFeatures; ++featureIdx)
    {
      importance[featureIdx] *= 1.0*nNodesInForest/ntrees_; //nContrastsInForest
      importance[featureIdx] /= nOobSamples; //nRealFeatures
    }

}

bool Randomforest::is_feature_in_tree(size_t featureidx, size_t treeidx)
{
  for(size_t i = 0; i < nnodes_[treeidx]; ++i)
    {
      if(forest_[treeidx][i].has_children())
	{
	  if(featureidx == forest_[treeidx][i].get_splitter())
	    {
	      return(true);
	    }
	}
    }
  return(false);
}

void Randomforest::tree_impurity(map<Node*,vector<size_t> >& trainics, num_t& impurity)
{

  impurity = 0.0;
  size_t n_tot = 0;
  
  size_t targetIdx = treedata_->getTarget();
  bool isTargetNumerical = treedata_->isFeatureNumerical(targetIdx);

  for(map<Node*,vector<size_t> >::iterator it(trainics.begin()); it != trainics.end(); ++it)
    {

      vector<num_t> targetData;
      treedata_->getFeatureData(targetIdx,it->second,targetData);

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
