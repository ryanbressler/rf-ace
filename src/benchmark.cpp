#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<vector>
#include<string>
#include "node.hpp"
#include "datadefs.hpp"

using namespace std;
using datadefs::cat_t;
using datadefs::num_t;

int main()
{
  int nsamples = 20;
  bool isregr = true;
  Node node_regr(nsamples,isregr);
  Node node_class(nsamples,!isregr);

  node_regr.print();
  node_class.print();

  int splitter = 6;
  num_t threshold = 3.4;
  int leftchild = 1;
  int rightchild = 2;
  node_regr.set_splitter(splitter,threshold,leftchild,rightchild);

  splitter = 11;
  set<cat_t> classet;
  classet.insert(1);
  classet.insert(2);
  classet.insert(4);
  leftchild = 3;
  rightchild = 5;
  node_class.set_splitter(splitter,classet,leftchild,rightchild);

  for(int i = 0; i < 15; ++i)
    {
      node_class.add_trainsample_idx(i);
      node_class.add_testsample_idx(2*i);
      node_regr.add_trainsample_idx(3*i);
      node_regr.add_testsample_idx(4*i);
    }

  
  node_regr.print();
  node_class.print();

  return(EXIT_SUCCESS);
}
