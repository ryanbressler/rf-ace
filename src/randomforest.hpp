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

  //Selects the target feature that is to be predicted
  void select_target(size_t targetidx);

  //Gets the selected target feature
  size_t get_target();

  //Grow the Random Forest with respect to selected target feature
  void grow_forest();

  //Grow the Random Forest with respect to targetidx
  void grow_forest(size_t targetidx);

private:

  //Grows one tree with respect to selected target feature
  void grow_tree(size_t treeidx);

  //Recursive tree-generating node splitter algorithm
  //NOTE: there will be at least two alternative node splitter algorithms in the future
  void recursive_nodesplit(size_t treeidx, size_t nodeidx, vector<size_t> const& sampleics);

  //Pointer to treedata_ object, stores all the feature information
  Treedata* treedata_;
  
  //Size of the forest
  size_t ntrees_;
  size_t mtry_;
  size_t nodesize_;

  vector<vector<Node> > forest_; //forest_[i][j] is the j'th node of i'th tree. forest_[i][0] is the rootnode.
  vector<size_t> nnodes_; //Number of used nodes in each tree.
  vector<vector<size_t> > oob_mat_;
  vector<size_t> noob_;

};

#endif
