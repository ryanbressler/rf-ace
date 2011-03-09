#ifndef RANDOMFOREST_HPP
#define RANDOMFOREST_HPP

#include<cstdlib>
#include "node.hpp"
#include "treedata.hpp"

using namespace std;

class Randomforest
{
public:

  //Grows a Random Forest with ntrees, mtry, and nodesize
  Randomforest(Treedata* treedata, size_t ntrees, size_t mtry, size_t nodesize);
  ~Randomforest();

private:

  void grow_forest();
  void grow_tree(size_t treeidx);
  void recursive_nodesplit(size_t treeidx, size_t nodeidx, vector<size_t> const& sampleics);

  Treedata* treedata_;
  
  size_t ntrees_;
  size_t mtry_;
  size_t nodesize_;

  vector<vector<Node> > forest_; //forest_[i][j] is the j'th node of i'th tree. forest_[i][0] is the rootnode.
  vector<size_t> nnodes_; //Number of used nodes in each tree.
  vector<vector<size_t> > oob_mat_;
  vector<size_t> noob_;

};

#endif
