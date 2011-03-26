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
//using datadefs::cat_t;
using datadefs::num_t;

class Node {
public:
  //Initializes node.
  Node();
  ~Node();

  //Sets a splitter feature for the node.
  //NOTE: splitter can be assigned only once! Subsequent setter calls will raise an assertion failure.
  void set_splitter(size_t splitter, set<num_t> classet, Node& leftchild, Node& rightchild);
  void set_splitter(size_t splitter, num_t threshold, Node& leftchild, Node& rightchild);

  //Gets the splitter for the node
  size_t get_splitter();
  
  //Given value, descends to either one of the child nodes if existent and returns true, otherwise false.
  //NOTE: childp is a ref-to-ptr that will be modified to point to the child node if descend is successful. 
  Node* percolate(num_t value);

  //THESE WILL POSSIBLY BECOME DEPRECATED
  void set_impurity(num_t value);
  num_t get_impurity();
  void reset_impurity();

  void set_trainidx(size_t trainidx);
  vector<size_t>* get_trainics();
  void clear_trainics();

  //Logic test whether the node has children or not
  bool has_children();

  //Helper functions
  void print();
  void print_compact();

private:
  bool isnum_;

  size_t splitter_;
  num_t threshold_;
  set<num_t> classet_;

  num_t impurity_;
    
  bool haschildren_;
  Node* leftchild_;
  Node* rightchild_;

  vector<size_t> trainics_;
};

#endif
