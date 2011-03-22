#include "randomforest.hpp"
#include<cmath>
#include<iostream>

Randomforest::Randomforest(Treedata* treedata, size_t ntrees, size_t mtry, size_t nodesize):
  treedata_(treedata),
  ntrees_(ntrees),
  mtry_(mtry),
  nodesize_(nodesize)
{

  size_t nsamples(treedata_->nsamples());

  //First we count the theoretical maximum number of nodes per tree.
  //Because each leaf must contain at least nodesize amount of data points, nmaxleaves is
  int nmaxleaves = int(ceil(float(nsamples)/nodesize_));
  //The upper bound for depth of the tree is log2(nmaxleaves)=log10(nmaxleaves)/log10(2.0):
  int maxdepth = int(ceil(log10(float(nmaxleaves))/log10(2.0)));
  //Thus, the number of nodes in a complete binary tree of depth maxdepth, there are
  int nmaxnodes = int(pow(2.0,maxdepth+1)); //In reality it's nmaxnodes-1 but this way we'll get a power of two which is supposed to be faster :)
  
  //Reserve memory for one tree
  vector<Node> tree(nmaxnodes);
 
  //Reserve memory for nnodes
  vector<size_t> nnodes(ntrees_);
  nnodes_ = nnodes;

  //Reserve memory for the whole forest
  vector<vector<Node> > forest(ntrees_);
  forest_ = forest;
  for(size_t i = 0; i < ntrees_; ++i)
    {
      forest_[i] = tree;
      
    }

  size_t defaulttargetidx = 0;
  Randomforest::select_target(defaulttargetidx);

  cout << forest_.size() << " trees and " << nmaxnodes << " max nodes per tree generated." << endl;

}

Randomforest::~Randomforest()
{

} 

void Randomforest::select_target(size_t targetidx)
{
  //cout << treedata_->get_target() << "\t" << targetidx << endl;
  if(treedata_->get_target() != targetidx)
    {
      treedata_->select_target(targetidx);
    }

  //Initialize oob_mat_ to store the out-of-box samples for the selected target
  size_t nrealvalues(treedata_->nrealvalues());
  vector<vector<size_t> > oob_mat(ntrees_);
  oob_mat_ = oob_mat;
  for(size_t i = 0; i < ntrees_; ++i)
    {
      vector<size_t> oob_ics(nrealvalues);
      oob_mat_[i] = oob_ics;
    }
  vector<size_t> noob(ntrees_);
  noob_ = noob;

  //cout << "Feature " << targetidx << " selected as target. Data sorted." << endl;
  //treedata_->print();
}

size_t Randomforest::get_target()
{
  return(treedata_->get_target());
}

void Randomforest::grow_forest()
{
  for(size_t i = 0; i < ntrees_; ++i)
    {
      Randomforest::grow_tree(i);
    }
}

void Randomforest::grow_forest(size_t targetidx)
{
  Randomforest::select_target(targetidx);
  Randomforest::grow_forest();
}

void Randomforest::grow_tree(size_t treeidx)
{
  //Generate the vector for bootstrap indices
  vector<size_t> bootstrap_ics(treedata_->nrealvalues());

  //Generate bootstrap indices, oob-indices, and noob
  treedata_->bootstrap(bootstrap_ics,oob_mat_[treeidx],noob_[treeidx]);

  //cout << "Growing tree " << treeidx << " with bootstrap sample:";
  //for(size_t i = 0; i < treedata_->nrealvalues(); ++i)
  //  {
  //    cout << " " << bootstrap_ics[i];
  //  }
  //cout << endl;

  size_t rootnode = 0;
  
  //Start the recursive node splitting from the root node. This will generate the tree.
  Randomforest::recursive_nodesplit(treeidx,rootnode,bootstrap_ics);
}

void Randomforest::recursive_nodesplit(size_t treeidx, size_t nodeidx, vector<size_t>& sampleics)
{

  size_t n_tot(sampleics.size());

  if(n_tot < 2*nodesize_)
    {
      return;
    }

  //Create mtry randomly selected feature indices to determine the split
  vector<size_t> mtrysample(treedata_->nfeatures());
  treedata_->permute(mtrysample);

  cout << "tree=" << treeidx << "  nodeidx=" << nodeidx << "  #sampleics=" << sampleics.size();
  
  vector<size_t> sampleics_left,sampleics_right;
  treedata_->split_target(nodesize_,sampleics,sampleics_left,sampleics_right);
  
  cout << "  #left=" << sampleics_left.size() << "  #right=" << sampleics_right.size() << endl;

  size_t n_left(sampleics_left.size());
  size_t n_right(sampleics_right.size());

  num_t bestrelativedecrease(0);
  size_t bestsplitter_i(mtry_);
  size_t targetidx(treedata_->get_target());
  for(size_t i = 0; i < mtry_; ++i)
    {

      size_t featureidx = mtrysample[i];

      if(featureidx == targetidx)
	{
	  continue;
	}
      
      num_t impurity_tot,impurity_left,impurity_right;

      treedata_->impurity(featureidx,sampleics,impurity_tot);
      treedata_->impurity(featureidx,sampleics_left,impurity_left);
      treedata_->impurity(featureidx,sampleics_right,impurity_right);

      num_t relativedecrease((impurity_tot-n_left*impurity_left/n_tot-n_right*impurity_right/n_tot)/impurity_tot);

      if(relativedecrease > bestrelativedecrease)
	{
	  bestrelativedecrease = relativedecrease;
	  bestsplitter_i = i;
	}
    }

  if(bestsplitter_i == mtry_)
    {
      cout << "No splitter found, quitting." << endl;
      return;
    }

  size_t splitterfeatureidx(mtrysample[bestsplitter_i]);

  cout << "Best splitter featureidx=" << splitterfeatureidx << " with relative decrease in impurity of " << bestrelativedecrease << endl; 
  
  //WE STILL HAVE TO MAKE THE FINAL SPLIT


  //treedata_->split_target_wrt_feature(splitterfeatureidx,nodesize_,sampleics,sampleics_left,sampleics_right);

  //featurematrix_[splitterfeatureidx]
  //featurematrix_[targetidx_]

  //vector<num_t> tv(n_tot);
  //vector<num_t> fv(n_tot);

  //for(size_t i = 0; i < n_tot; ++i)
  //  {
  //    tv[i] = featurematrix_[targetidx_][sampleics[i]];
  //    fv[i] = featurematrix_[splitterfeatureidx][sampleics[i]];
  //  }

  //vector<size_t> ref_ics(n_tot);

  //These functions will be moved under data definitions
  //treedata_->sort_and_make_ref<num_t>(fv,ref_ics);
  //treedata_->sort_from_ref<num_t>(tv,ref_ics);
  //treedata_->sort_from_ref<size_t>(sampleics);

  size_t nodeidx_left(++nnodes_[treeidx]);
  size_t nodeidx_right(++nnodes_[treeidx]);

  forest_[treeidx][nodeidx].set_splitter(mtrysample[bestsplitter_i],datadefs::num_nan,forest_[treeidx][nodeidx_left],forest_[treeidx][nodeidx_right]);
  
  Randomforest::recursive_nodesplit(treeidx,nodeidx_left,sampleics_left);
  Randomforest::recursive_nodesplit(treeidx,nodeidx_right,sampleics_right);
  
  
}
