#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<vector>
#include<string>
#include "node.hpp"
#include "datadefs.hpp"
#include "treedata.hpp"

using namespace std;
using datadefs::cat_t;
using datadefs::num_t;

//This function will demonstrate the functionality of the current development stage of the rf-ace program. The function includes:
//DONE:generate node hierarchies using the node class (DONE)
//DONE:demonstrate percolation of samples through the tree using the descend() function (DONE)
//DONE:demonstrate adding sample indices to the nodes (DONE)
//NOT DONE:read mixed-type data with missing values into treedata object
//NOT DONE:embed artificial contrasts into treedata
//MORE TASKS TO COME 
int main()
{

  cout << endl;
  cout << "---------------------------------------" << endl;
  cout << "PART 1: DEMONSTRATE USAGE OF NODE CLASS" << endl;
  cout << "-create nodes" << endl;
  cout << "-add hierarchies and splitters" << endl;
  cout << "-percolate samples" << endl;
  cout << "---------------------------------------" << endl << endl;
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
      if(rootnode_cat.has_children())
	{
	  childp = rootnode_cat.percolate(i);
	  childp->add_trainsample_idx(i);
	}
      num_t j = i*1.0;
      if(rootnode_num.has_children())
      	{
	  childp = rootnode_num.percolate(j);
	  childp->add_trainsample_idx(i);
	}
   }

  rootnode_num.print();
  leftchild_num.print();
  rightchild_num.print();

  rootnode_cat.print();
  leftchild_cat.print();
  rightchild_cat.print();

  cout << endl;
  cout << "-------------------------------------------" << endl;
  cout << "PART 2: DEMONSTRATE USAGE OF TREEDATA CLASS" << endl;
  cout << "-..." << endl;
  cout << "-..." << endl;
  cout << "-..." << endl;
  cout << "-------------------------------------------" << endl << endl;



  return(EXIT_SUCCESS);
}
