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
  void setSplitter(size_t splitter, set<num_t> classSet, Node& leftChild, Node& rightChild);
  void setSplitter(size_t splitter, num_t threshold, Node& leftChild, Node& rightChild);

  //Gets the splitter for the node
  inline size_t getSplitter() { return(splitter_); }

  //Given value, descends to either one of the child nodes if existent and returns true, otherwise false.
  //NOTE: childp is a ref-to-ptr that will be modified to point to the child node if descend is successful. 
  Node* percolateData(num_t value);

  //THESE WILL POSSIBLY BECOME DEPRECATED
  //void set_impurity(num_t value);
  //num_t get_impurity();
  //void reset_impurity();

  //void set_prediction(num_t value);
  //num_t get_prediction();

  //void set_trainidx(size_t trainidx);
  //vector<size_t>* get_trainics();
  //void clearTrainIcs();

  //Logic test whether the node has children or not
  inline bool hasChildren() { return(hasChildren_); }

  //Helper functions
  void print();
  void print_compact();

private:
  bool isSplitterNumerical_;

  size_t splitter_;
  num_t threshold_;
  set<num_t> classSet_;

  //num_t impurity_;
  //num_t prediction_; // saves the prediction of the node
    
  bool hasChildren_;
  Node* leftChild_;
  Node* rightChild_;

  //vector<size_t> trainIcs_;
};

#endif
