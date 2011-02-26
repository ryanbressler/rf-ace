//node.hpp
//
//A node class for CARTs

#ifndef NODE_HPP
#define NODE_HPP

#include<cstdlib>
#include<vector>
#include<set>
#include "datadefs.hpp"

using namespace std;
using datadefs::cat_t;
using datadefs::num_t;

class Node {
public:
  //Initializes node to store at max nsamples, either numerical (isregr == true) or categorical (isregr == false).
  //Excess memory will be reserved in order to avoid dynamic memory allocation.
  Node(int nsamples, bool isregr);
  ~Node();

  //Sets a splitter feature for the node.
  //NOTE: splitter can be assigned only once! Subsequent setter calls will raise an assertion failure.
  void set_splitter(int splitter, set<cat_t> classet, int leftchild, int rightchild);
  void set_splitter(int splitter, num_t threshold, int leftchild, int rightchild);

  //Gets the splitter feature for the node
  int get_splitter();
  
  //Given value, descends to either one of the child nodes if existent and returns true, otherwise false.
  //NOTE: the nodes aren't aware of each other, they just know the indices referring to their children, if existent.
  //NOTE: thus, nodes don't constitute a tree per se, that must be performed inside the CART implementation.
  bool descend(cat_t value, int& child);
  bool descend(num_t value, int& child);

  void add_trainsample_idx(int idx);
  void add_testsample_idx(int idx);

  void reset_testsample_ics();

  //Logic test whether the node has children or not
  bool is_leaf();

  //Helper functions
  void print();
  void print_compact();

private:
  bool isregr_;

  int splitter_;
  num_t threshold_;
  set<cat_t> classet_;

  vector<int> trainsampleics_;
  vector<int> testsampleics_;

  //vector<num_t> num_trainsamples_;
  //vector<num_t> num_testsamples_;

  size_t ntrainsamples_;
  size_t ntestsamples_;
  
  bool haschildren_;
  int leftchild_;
  int rightchild_;
};

#endif
