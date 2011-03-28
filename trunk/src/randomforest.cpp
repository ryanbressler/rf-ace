#include "randomforest.hpp"
#include "datadefs.hpp"
#include<cmath>
#include<iostream>

Randomforest::Randomforest(Treedata* treedata, size_t ntrees, size_t mtry, size_t nodesize):
  treedata_(treedata),
  ntrees_(ntrees),
  mtry_(mtry),
  nodesize_(nodesize),
  oobmatrix_(ntrees),
  trainics_(ntrees)
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

  //for(size_t i = 0; i < ntrees_; ++i)
  // {
  //   vector<size_t> foo;
  //   oobmatrix_[i] = foo;
  // }

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

  for(size_t i = 0; i < ntrees_; ++i)
    {
      oobmatrix_[i].clear();
    }
  //for(size_t i = 0; i < ntrees_; ++i)
  //  {
  //    vector<size_t> oob_ics(nrealvalues);
  //    oob_mat_[i] = oob_ics;
  //  }
  //vector<size_t> noob(ntrees_);
  //noob_ = noob;

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
      cout << "starting to grow treeidx=" << i << endl;
      Randomforest::grow_tree(i);
    }
}

//void Randomforest::grow_forest(size_t targetidx)
//{
//  Randomforest::select_target(targetidx);
//  Randomforest::grow_forest();
//}

void Randomforest::grow_tree(size_t treeidx)
{
  //Generate the vector for bootstrap indices
  vector<size_t> bootstrap_ics(treedata_->nrealvalues());

  //Generate bootstrap indices and oob-indices
  treedata_->bootstrap(bootstrap_ics,oobmatrix_[treeidx]);

  //THIS IS FOR TESTING PURPOSES ONLY
  //vector<size_t> bootstrap_ics(5);
  //vector<size_t> foo(5);
  //for(size_t i = 0; i < 5; ++i)
  //  {
  //    foo[i] = i+5;
  //    bootstrap_ics[i] = i;
  //  }
  //cout << "right before copying to oobmatrix_...";
  //oobmatrix_[treeidx] = foo;
  //cout << " done" << endl;

  cout << "tree=" << treeidx << "  bootstrap indices [";
  for(size_t i = 0; i < bootstrap_ics.size(); ++i)
    {
      cout << " " << bootstrap_ics[i];
    }
  cout << " ]  oob [";
  for(size_t i = 0; i < oobmatrix_[treeidx].size(); ++i)
    {
      cout << " " << oobmatrix_[treeidx][i];
    }
  cout << " ]" << endl;

  size_t rootnode = 0;
  
  //Start the recursive node splitting from the root node. This will generate the tree.
  Randomforest::recursive_nodesplit(treeidx,rootnode,bootstrap_ics);
  
  //Percolate oob samples
  Randomforest::percolate_sampleics(treeidx,oobmatrix_[treeidx]);
}

void Randomforest::recursive_nodesplit(size_t treeidx, size_t nodeidx, vector<size_t>& sampleics)
{

  //return;

  size_t n_tot(sampleics.size());

  if(n_tot < 2*nodesize_)
    {
      return;
    }

  //Create mtry randomly selected feature indices to determine the split
  vector<size_t> mtrysample(treedata_->nfeatures());
  treedata_->permute(mtrysample);

  cout << "tree=" << treeidx << "  nodeidx=" << nodeidx;

  //This is testing...
  vector<size_t> sampleics_left,sampleics_right;
  set<num_t> values_left;
  ///HERE'S SOMETHING GOING WRONG
  //treedata_->split_target_with_cat_feature(4,nodesize_,sampleics,sampleics_left,sampleics_right,values_left);
  //cout << "jee" << endl;
  treedata_->split_target(nodesize_,sampleics,sampleics_left,sampleics_right);

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
      
      if(impurity_tot < datadefs::eps)
	{
	  continue;
	}

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

  //cout << "Best splitter featureidx=" << splitterfeatureidx << " with relative decrease in impurity of " << bestrelativedecrease << endl; 

  size_t nodeidx_left(++nnodes_[treeidx]);
  size_t nodeidx_right(++nnodes_[treeidx]);

  if(treedata_->isfeaturenum(splitterfeatureidx))
    {
      num_t splitvalue;
      treedata_->split_target_with_num_feature(splitterfeatureidx,nodesize_,sampleics,sampleics_left,sampleics_right,splitvalue);
      if(sampleics_left.size() == 0 || sampleics_right.size() == 0)
	{
	  return;
	}
      forest_[treeidx][nodeidx].set_splitter(splitterfeatureidx,splitvalue,forest_[treeidx][nodeidx_left],forest_[treeidx][nodeidx_right]);
      //cout << forest_[treeidx][nodeidx].has_children() << endl;
    }
  else
    {
      set<num_t> values_left;
      treedata_->split_target_with_cat_feature(splitterfeatureidx,nodesize_,sampleics,sampleics_left,sampleics_right,values_left);
      if(sampleics_left.size() == 0 || sampleics_right.size() == 0)
	{
	  return;
	}
      forest_[treeidx][nodeidx].set_splitter(splitterfeatureidx,values_left,forest_[treeidx][nodeidx_left],forest_[treeidx][nodeidx_right]);
      //cout << forest_[treeidx][nodeidx].has_children() << endl;
    }
  
  Randomforest::recursive_nodesplit(treeidx,nodeidx_left,sampleics_left);
  Randomforest::recursive_nodesplit(treeidx,nodeidx_right,sampleics_right);
  
  
}


void Randomforest::percolate_sampleics(size_t treeidx, vector<size_t>& sampleics)
{
  
  for(size_t i = 0; i < sampleics.size(); ++i)
    {
      Node* nodep(&forest_[treeidx][0]);
      size_t sampleidx(sampleics[i]);
      treedata_->percolate_sampleidx(sampleidx,&nodep);
      map<Node*,vector<size_t> >::iterator it(trainics_[treeidx].find(nodep));
      if(it == trainics_[treeidx].end())
	{
	  Node* foop(nodep);
	  vector<size_t> foo(1);
	  foo[0] = sampleidx;
	  trainics_[treeidx].insert(pair<Node*,vector<size_t> >(foop,foo));
	}
      else
      	{
      	  trainics_[treeidx][it->first].push_back(sampleidx);
      	}
      
    }
  
  cout << "treeidx=" << treeidx << ", train samples percolated accordingly:" << endl;
  size_t iter = 0;
  for(map<Node*,vector<size_t> >::const_iterator it(trainics_[treeidx].begin()); it != trainics_[treeidx].end(); ++it, ++iter)
    {
      cout << "leaf node " << iter << ":"; 
      for(size_t i = 0; i < it->second.size(); ++i)
	{
	  cout << " " << it->second[i];
      	}
      cout << endl;
    }
  //cout << "done" << endl;
  
}

