#include "randomforest.hpp"
#include "datadefs.hpp"
#include<cmath>
#include<iostream>

Randomforest::Randomforest(Treedata* treedata, size_t ntrees, size_t mtry, size_t nodesize):
  treedata_(treedata),
  ntrees_(ntrees),
  mtry_(mtry),
  nodesize_(nodesize),
  oobmatrix_(ntrees)
{

  size_t defaulttargetidx = 0;
  Randomforest::select_target(defaulttargetidx);

  cout << forest_.size() << " trees and " << forest_[0].size() << " max nodes per tree initialized." << endl;

}

Randomforest::~Randomforest()
{

} 

void Randomforest::init_forest()
{

  size_t nsamples = treedata_->nsamples();

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
    }
}
  
void Randomforest::select_target(size_t targetidx)
{
  if(treedata_->get_target() != targetidx)
    {
      treedata_->select_target(targetidx);
    }

  for(size_t i = 0; i < ntrees_; ++i)
    {
      oobmatrix_[i].clear();
    } 
  Randomforest::init_forest();
}

size_t Randomforest::get_target()
{
  return(treedata_->get_target());
}

void Randomforest::grow_forest(const size_t nperms, const num_t alpha, vector<num_t>& pvalues)
{
  assert(nperms > 5);
  vector<vector<num_t> > importancemat(nperms);
  //vector<num_t> csample(nperms);

  for(size_t p = 0; p < nperms; ++p)
    {
      cout << "Growing forest " << p << endl;
      Randomforest::init_forest();
      treedata_->permute_contrasts();
      for(size_t t = 0; t < ntrees_; ++t)
	{
	  Randomforest::grow_tree(t);
	}
      //Randomforest::calculate_importance(alpha,importancemat[p],csample[p]);
      Randomforest::calculate_importance(importancemat[p]);
    }

  size_t nfeatures = treedata_->nfeatures();
  pvalues.resize(nfeatures);
  for(size_t f = 0; f < nfeatures; ++f)
    {
      vector<num_t> fsample(nperms);
      vector<num_t> csample(nperms);
      for(size_t p = 0; p < nperms; ++p)
	{
	  fsample[p] = importancemat[p][f];
	  csample[p] = importancemat[p][f + nfeatures];
	}
      datadefs::ttest(fsample,csample,pvalues[f]);
    }
}

void Randomforest::grow_tree(size_t treeidx)
{
  //Generate the vector for bootstrap indices
  vector<size_t> bootstrap_ics(treedata_->nrealvalues());

  //Generate bootstrap indices and oob-indices
  treedata_->bootstrap(bootstrap_ics,oobmatrix_[treeidx]);

  if(false)
    {
      cout << "tree " << treeidx << "  bootstrap indices [";
      for(size_t i = 0; i < bootstrap_ics.size(); ++i)
	{
	  cout << " " << bootstrap_ics[i];
	}
      cout << " ]  oob [";
      for(size_t i = 0; i < oobmatrix_[treeidx].size(); ++i)
	{
	  cout << " " << oobmatrix_[treeidx][i];
	}
      cout << " ]" << endl << endl;
    }

  size_t rootnode = 0;
  
  //Start the recursive node splitting from the root node. This will generate the tree.
  Randomforest::recursive_nodesplit(treeidx,rootnode,bootstrap_ics);
  
  cout << "Tree " << treeidx << ", nodes:";
  for(size_t i = 0; i < nnodes_[treeidx]; ++i)
    {
      cout << "|";
    }
  cout << endl;

}

void Randomforest::recursive_nodesplit(size_t treeidx, size_t nodeidx, vector<size_t>& sampleics)
{

  size_t n_tot(sampleics.size());

  if(n_tot < 2*nodesize_)
    {
      //cout << "Too few samples to start with, quitting" << endl;
      return;
    }

  //size_t nfeatures = treedata_->nfeatures();
  //Create mtry randomly selected feature indices to determine the split
  vector<size_t> mtrysample(mtry_);

  //size_t contrast_lower_limit = mtry_/2;
  
  for(size_t i = 0; i < mtry_; ++i)
    {
      size_t nallfeatures = 2*treedata_->nfeatures();
      mtrysample[i] = treedata_->randidx(nallfeatures);
    }

  //vector<bool> iscontrast(mtry_);
  //treedata_->permute(mtrysample);
  //treedata_->generate_contrasts(iscontrast);

  //cout << "Tree " << treeidx << "  Node " << nodeidx << endl;

  vector<size_t> sampleics_left,sampleics_right;
  num_t splitvalue;
  set<num_t> values_left;
  treedata_->split_target(treedata_->get_target(),nodesize_,sampleics,sampleics_left,sampleics_right,splitvalue,values_left);

  assert(n_tot == sampleics.size());

  //size_t nreal_tot;

  //num_t fitness = 0.0;
  num_t bestfitness = 0.0;
  size_t bestfeatureidx = mtry_;
  size_t targetidx = treedata_->get_target();

  //vector<num_t> fitness(mtry_);

  for(size_t i = 0; i < mtry_; ++i)
    {
      size_t featureidx = mtrysample[i];

      if(featureidx == targetidx)
        {
          continue;
        }

      num_t fitness = treedata_->split_fitness(featureidx,nodesize_,sampleics,sampleics_left,sampleics_right);

      if(fitness > bestfitness)
	{
	  bestfitness = fitness;
	  bestfeatureidx = featureidx;
	  if(fabs(bestfitness - 1.0) < datadefs::eps)
	    {
	      //cout << "Maximum fitness reached for splitter " << bestfeatureidx << ", stopped searching" << endl;
	      break;
	    }
	}
    }
  
  if(bestfeatureidx == mtry_)
    {
      //cout << "No splitter found, quitting" << endl << endl;
      return;
    }
  
  //cout << "Best splitter feature is " << splitterfeatureidx << " with relative decrease in impurity of " << bestrelativedecrease << endl; 
  
  //vector<size_t> sampleics_copy = sampleics;

  size_t nreal_tot;
  treedata_->remove_nans(bestfeatureidx,sampleics,nreal_tot);
  //cout << "Splitter " << bestfeatureidx << " has " << n_tot - nreal_tot << " missing values, which will be omitted in splitting" << endl;
  n_tot = nreal_tot;

  if(n_tot < 2*nodesize_)
    {
      //cout << "Splitter has too few non-missing values, quitting" << endl;
      //This needs to be fixed such that one of the surrogates will determine the split instead
      return;
    }

  size_t nodeidx_left; //(++nnodes_[treeidx]);
  size_t nodeidx_right; //(++nnodes_[treeidx]);

  //if(treedata_->isfeaturenum(bestfeatureidx))
    // {
  //num_t splitvalue;
  //set<num_t> values_left;
  treedata_->split_target(bestfeatureidx,nodesize_,sampleics,sampleics_left,sampleics_right,splitvalue,values_left);
  assert(sampleics.size() == n_tot);
  assert(sampleics_left.size() + sampleics_right.size() == n_tot);
  if(sampleics_left.size() < nodesize_ || sampleics_right.size() < nodesize_)
    {
      //cout << "Split was unsuccessful, quitting" << endl;
      return;
    }
  
  nodeidx_left = ++nnodes_[treeidx];
  nodeidx_right = ++nnodes_[treeidx];
  
  if(treedata_->isfeaturenum(bestfeatureidx))
    {
      forest_[treeidx][nodeidx].set_splitter(bestfeatureidx,splitvalue,forest_[treeidx][nodeidx_left],forest_[treeidx][nodeidx_right]);
    }
  else
    {
      forest_[treeidx][nodeidx].set_splitter(bestfeatureidx,values_left,forest_[treeidx][nodeidx_left],forest_[treeidx][nodeidx_right]);
    }
  
  if(sampleics_left.size() > 2*nodesize_)
    {
      Randomforest::recursive_nodesplit(treeidx,nodeidx_left,sampleics_left);
    }

  if(sampleics_right.size() > 2*nodesize_)
    {
      Randomforest::recursive_nodesplit(treeidx,nodeidx_right,sampleics_right);
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
  
  
  /*
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
  */
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
      num_t value(treedata_->at(featureidx_new,sampleidx));
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
	  value = treedata_->randf(featureidx_new);
	}
      else
	{
	  value = treedata_->at(featureidx_new,sampleidx);
	}
      *nodep = (*nodep)->percolate(value);
    }
}


void Randomforest::calculate_importance(vector<num_t>& importance)
{

  size_t nrealfeatures(treedata_->nfeatures());
  size_t nallfeatures(2*nrealfeatures);
  importance.resize(nallfeatures);
  size_t noob_tot(0);

  for(size_t i = 0; i < nallfeatures; ++i)
    {
      importance[i] = 0;
    }

  for(size_t i = 0; i < ntrees_; ++i)
    {
      size_t noob_new(oobmatrix_[i].size());
      noob_tot += noob_new;
      Node rootnode(forest_[i][0]);
      map<Node*,vector<size_t> > trainics;
      Randomforest::percolate_sampleics(rootnode,oobmatrix_[i],trainics);
      num_t impurity_tree;
      Randomforest::tree_impurity(trainics,impurity_tree);
      //cout << "#nodes_with_train_samples=" << trainics.size() << endl;  
      for(size_t f = 0; f < nallfeatures; ++f)
	{
	  if(Randomforest::is_feature_in_tree(f,i))
	    {
	      Randomforest::percolate_sampleics_randf(f,rootnode,oobmatrix_[i],trainics);
	      num_t impurity_perm;
	      Randomforest::tree_impurity(trainics,impurity_perm);
	      importance[f] += noob_new * (impurity_perm - impurity_tree) / impurity_tree; 
	    }
	}
    }

  
  for(size_t i = 0; i < nallfeatures; ++i)
    {
      importance[i] /= noob_tot;
    }
    
  /*
    vector<num_t> contrast_importance(nrealfeatures);
    for(size_t i = nrealfeatures; i < nallfeatures; ++i)
    {
    contrast_importance[i - nrealfeatures] = importance[i];
    }
    importance.resize(nrealfeatures);
    
    datadefs::percentile(contrast_importance,alpha,contrast_prc);
    cout << "CONTRAST: " << 100*alpha << "th percentile is " << contrast_prc << endl;
    
    if(false)
    {
    for(size_t i = 0; i < nrealfeatures; ++i)
    {
    if(importance[i] > contrast_prc)
    {
    cout << treedata_->get_targetheader() << "\t" << treedata_->get_featureheader(i) << "\t" << importance[i] << endl;
    }
    }
    }
  */
  

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
  size_t n_tot(0);
  
  size_t targetidx(treedata_->get_target());

  for(map<Node*,vector<size_t> >::iterator it(trainics.begin()); it != trainics.end(); ++it)
    {
      num_t impurity_leaf;
      size_t nreal;
      treedata_->impurity(targetidx,it->second,impurity_leaf,nreal);
      //size_t n(it->second.size());
      n_tot += nreal;
      impurity += nreal * impurity_leaf;
    }

  if(n_tot > 0)
    {
      impurity /= n_tot;
    }
}
