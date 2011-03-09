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

  //Reserve memory for the oob-ics and noob
  vector<vector<size_t> > oob_mat(ntrees_);
  vector<size_t> oob_ics(nsamples);
  for(size_t i = 0; i < ntrees_; ++i)
    {
      oob_mat[i] = oob_ics;
    }
  oob_mat_ = oob_mat;
  vector<size_t> noob(ntrees_);
  noob_ = noob;

  cout << "Forest initialized. " << forest_.size() << " trees and " 
       << nmaxnodes << " max nodes per tree generated." << endl;

  //Finally, grow the forest
  Randomforest::grow_forest();

}

Randomforest::~Randomforest()
{

}

void Randomforest::grow_forest()
{
  for(size_t i = 0; i < ntrees_; ++i)
    {
      Randomforest::grow_tree(i);
    }
}

void Randomforest::grow_tree(size_t treeidx)
{
  size_t nsamples(treedata_->nsamples());
  //Generate the vector for bootstrap (sample) indices
  vector<size_t> bootstrap_ics(nsamples);

  //Generate bootstrap indices, oob-indices, and noob
  treedata_->bootstrap(bootstrap_ics,oob_mat_[treeidx],noob_[treeidx]);

  cout << "Growing tree " << treeidx << " with bootstrap sample:";
  for(size_t i = 0; i < nsamples; ++i)
    {
      cout << " " << bootstrap_ics[i];
    }
  cout << endl;

  size_t rootnode = 0;
  
  //Start the recursive node splitting from the root node. This will generate the tree.
  Randomforest::recursive_nodesplit(treeidx,rootnode,bootstrap_ics);
}

void Randomforest::recursive_nodesplit(size_t treeidx, size_t nodeidx, vector<size_t> const& sampleics)
{

  cout << "Splitting..." << endl;
  forest_[treeidx][nodeidx].print();

}
