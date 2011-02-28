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

  Node rootnode_num(nsamples);
  Node leftchild_num(nsamples);
  Node rightchild_num(nsamples);

  Node rootnode_cat(nsamples);
  Node leftchild_cat(nsamples);
  Node rightchild_cat(nsamples);

  int splitter_num = 6;
  num_t threshold = 3.4;
  rootnode_num.set_splitter(splitter_num,threshold,leftchild_num,rightchild_num);

  int splitter_cat = 11;
  set<cat_t> classet;
  classet.insert(1);
  classet.insert(2);
  classet.insert(4);
  rootnode_cat.set_splitter(splitter_cat,classet,leftchild_cat,rightchild_cat);

  Node* childp(NULL); 
  for(int i = 0; i < 10; ++i)
    {
      if(rootnode_cat.descend(i,&childp))
	{
	  childp->add_trainsample_idx(i);
	}
      num_t j = i*1.0;
      if(rootnode_num.descend(j,&childp))
      	{
	  childp->add_trainsample_idx(i);
	}
   }

  rootnode_num.print();
  leftchild_num.print();
  rightchild_num.print();

  rootnode_cat.print();
  leftchild_cat.print();
  rightchild_cat.print();

  return(EXIT_SUCCESS);
}
