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
//DONE:demonstrate percolation of samples through the tree (DONE)
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
  //cout << "-percolate samples" << endl;
  cout << "---------------------------------------" << endl << endl;
  //int nsamples = 20;

  Node node; //rootnode 
  Node node1; //left (LEAF)
  Node node2; //right
  Node node21; //right+left
  Node node22; //right+right (LEAF)
  Node node211; //right+left+left (LEAF)
  Node node212; //right+left+right (LEAF)

  int splitter_num = 6;
  num_t threshold = 3.4;
  node.set_splitter(splitter_num,threshold,node1,node2);

  int splitter_cat = 11;
  set<cat_t> classet;
  classet.insert(1);
  classet.insert(2);
  classet.insert(4);
  node2.set_splitter(splitter_cat,classet,node21,node22);
  
  splitter_num = 7;
  threshold = 5.6;
  node21.set_splitter(splitter_num,threshold,node211,node212);

  node.print();
  node1.print();
  node2.print();
  node21.print();
  node22.print();
  node211.print();
  node212.print();

  cout << endl;
  cout << "-------------------------------------------" << endl;
  cout << "PART 2: DEMONSTRATE USAGE OF TREEDATA CLASS" << endl;
  cout << "-Load mixed-type data with missing values" << endl;
  cout << "-..." << endl;
  cout << "-..." << endl;
  cout << "-------------------------------------------" << endl << endl;

  bool is_featurerows = true;
  string fname = "data/test_6by10_featurerows_matrix.tsv";
  Treedata treedata(fname,is_featurerows);
  treedata.select_target(0);
  treedata.print();

  return(EXIT_SUCCESS);
}
