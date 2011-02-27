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

  Node rootnode_regr(nsamples,isregr);
  Node leftchild_regr(nsamples,isregr);
  Node rightchild_regr(nsamples,isregr);

  Node rootnode_class(nsamples,!isregr);
  Node leftchild_class(nsamples,!isregr);
  Node rightchild_class(nsamples,!isregr);

  int splitter_regr = 6;
  num_t threshold = 3.4;
  rootnode_regr.set_splitter(splitter_regr,threshold,leftchild_regr,rightchild_regr);

  int splitter_class = 11;
  set<cat_t> classet;
  classet.insert(1);
  classet.insert(2);
  classet.insert(4);
  rootnode_class.set_splitter(splitter_class,classet,leftchild_class,rightchild_class);

  Node* childp(NULL); 
  for(int i = 0; i < 10; ++i)
    {
      if(rootnode_class.descend(i,&childp))
	{
	  childp->add_trainsample_idx(i);
	}
      num_t j = i*1.0;
      if(rootnode_regr.descend(j,&childp))
      	{
	  childp->add_trainsample_idx(i);
	}
   }

  rootnode_regr.print();
  leftchild_regr.print();
  rightchild_regr.print();

  rootnode_class.print();
  leftchild_class.print();
  rightchild_class.print();

  return(EXIT_SUCCESS);
}
