#ifndef RANDOMFOREST_HPP
#define RANDOMFOREST_HPP

#include<cstdlib>
#include "node.hpp"
#include "treedata.hpp"

using namespace std;

class Randomforest
{
public:
  Randomforest(Treedata* treedata, size_t ntrees, size_t mtry, size_t nodesize);
  ~Randomforest();

private:

  Treedata* treedata_;

  size_t ntrees_;
  size_t mtry_;
  size_t nodesize_;

};

#endif
