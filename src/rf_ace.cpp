#include<cstdlib>
#include "randomforest.hpp"
#include "treedata.hpp"

using namespace std;

int main()
{
  //FIRST PART: read data into Treedata class (features are rows)
  bool is_featurerows = true;
  string fname = "data/test_6by10_featurerows_matrix.tsv";
  Treedata treedata(fname,is_featurerows);
  
  //SECOND PART: construct a Random Forest object
  size_t ntrees(5);
  size_t mtry(3);
  size_t nodesize(1);
  Randomforest RF(&treedata,ntrees,mtry,nodesize);
  
  size_t targetidx(2);
  RF.grow_forest(targetidx);



  return(EXIT_SUCCESS);
}
