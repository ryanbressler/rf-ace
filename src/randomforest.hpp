#ifndef RANDOMFOREST_HPP
#define RANDOMFOREST_HPP

#include<cstdlib>
#include "node.hpp"
#include "treedata.hpp"

using namespace std;

class Randomforest
{
public:

  //Initializes Random Forest with ntrees, mtry, and nodesize
  //NOTE: target feature will be set to 0. Change with Randomforest::select_target()
  Randomforest(Treedata* treedata, size_t ntrees, size_t mtry, size_t nodesize);
  ~Randomforest();

  void init_forest();

  //Selects the target feature that is to be predicted
  void select_target(size_t targetidx);

  //Gets the selected target feature
  size_t get_target();

  //Grow the Random Forest with respect to selected target feature
  void grow_forest(const size_t nperms, vector<num_t>& pvalues, vector<num_t>& ivalues);
  //void calculate_importance(const num_t alpha, vector<num_t>& importance, num_t& contrast_prc);

  //blacklist_and_kill

private:

  //Grows one tree with respect to selected target feature
  void grow_tree(size_t treeidx);

  //Recursive tree-generating node splitter algorithm
  //NOTE: there will be at least two alternative node splitter algorithms in the future
  void recursive_nodesplit(size_t treeidx, size_t nodeidx, vector<size_t>& sampleics);

  void percolate_sampleics(Node& rootnode, vector<size_t>& sampleics, map<Node*,vector<size_t> >& trainics);
  void percolate_sampleics_randf(size_t featureidx, Node& rootnode, vector<size_t>& sampleics, map<Node*,vector<size_t> >& trainics);
  
  void percolate_sampleidx(size_t sampleidx, Node** nodep);
  void percolate_sampleidx_randf(size_t featureidx, size_t sampleidx, Node** nodep);

  bool is_feature_in_tree(size_t featureidx, size_t treeidx);

  void tree_impurity(map<Node*,vector<size_t> >& trainics, num_t& impurity);

  //void calculate_importance(const num_t alpha, vector<num_t>& importance, num_t& contrast_prc);
  void calculate_importance(vector<num_t>& importance);

  //Pointer to treedata_ object, stores all the feature information
  Treedata* treedata_;
  
  //Size of the forest
  size_t ntrees_;
  size_t mtry_;
  size_t nodesize_;

  vector<vector<Node> > forest_; //forest_[i][j] is the j'th node of i'th tree. forest_[i][0] is the rootnode.
  vector<size_t> nnodes_; //Number of used nodes in each tree.
  vector<vector<size_t> > oobmatrix_;
  
  //vector<map<Node*,vector<size_t> > > trainics_;

};

#endif
