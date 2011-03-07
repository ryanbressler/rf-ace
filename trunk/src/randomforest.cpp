#include "randomforest.hpp"

Randomforest::Randomforest(Treedata* treedata, size_t ntrees, size_t mtry, size_t nodesize):
  treedata_(treedata),
  ntrees_(ntrees),
  mtry_(mtry),
  nodesize_(nodesize)
{
  
}

Randomforest::~Randomforest()
{
  
}
