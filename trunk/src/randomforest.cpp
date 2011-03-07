#include "randomforest.hpp"
#include<cmath>
#include<iostream>

Randomforest::Randomforest(Treedata* treedata, size_t ntrees, size_t mtry, size_t nodesize):
  treedata_(treedata),
  ntrees_(ntrees),
  mtry_(mtry),
  nodesize_(nodesize)
{

  //First we count the theoretical maximum number of nodes per tree.
  //Because each leaf must contain at least nodesize amount of data points, nmaxleaves is
  int nmaxleaves = int(ceil(float(treedata_->nsamples())/nodesize_));
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

  cout << "Forest initialized. " << forest_.size() << " trees and " 
       << nmaxnodes << " max nodes per tree generated." << endl;

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
  forest_[treeidx][0].print();
}
