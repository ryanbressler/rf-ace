#include <cstdlib>
#include <cassert>
#include <string>
#include <iostream>
#include "randomforest.hpp"
#include "treedata.hpp"

using namespace std;

int main()
{
  
  Treedata treedata("data/big_real_matrix.tsv",true);
  Randomforest RF(&treedata,2,10,10);
  RF.select_target(0);
  RF.grow_forest();

  return(EXIT_SUCCESS);
}
