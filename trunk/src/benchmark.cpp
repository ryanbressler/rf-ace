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
  int nsamples = 5;
  bool isregr = true;
  Node node_regr(nsamples,isregr);
  Node node_class(nsamples,!isregr);

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

  cout << "Classification node: splitter is " << node_class.get_splitter() << endl;
  cout << "Regression node: splitter is " << node_regr.get_splitter() << endl;
  

  return(EXIT_SUCCESS);
}
