//node.hpp
//
//A node class for CARTs

#ifndef NODE_HPP
#define NODE_HPP

#include <cstdlib>
#include <vector>
#include <set>
#include "datadefs.hpp"

using namespace std;
using datadefs::num_t;

class Node {
public:
  //Initializes node.
  Node();
  ~Node();

  //Sets a splitter feature for the node.
  //NOTE: splitter can be assigned only once! Subsequent setter calls will raise an assertion failure.
  void setSplitter(size_t splitterIdx, set<num_t> classSet);
  void setSplitter(size_t splitterIdx, num_t threshold);

  //Gets the splitter for the node
  inline size_t getSplitter() { return(splitter_); }

  //Given a value, descends to either one of the child nodes, if existing, otherwise returns a pointer to the current node
  Node* percolateData(num_t value);

  void setPrediction(num_t value);
  num_t getPrediction();

  //void set_trainidx(size_t trainidx);
  //vector<size_t>* get_trainics();
  //void clearTrainIcs();

  //Logic test whether the node has children or not
  inline bool hasChildren() { return(hasChildren_); }

  Node* leftChild();
  Node* rightChild();

  size_t nNodes();

  //Helper functions
  //void print();
  //void print_compact();

private:

  //A recursive function that will accumulate the number of descendant nodes the current nodes has
  void recursiveNDescendantNodes(size_t& n);

  bool isSplitterNumerical_;

  size_t splitter_;
  num_t threshold_;
  set<num_t> classSet_;

  //num_t impurity_;
  num_t prediction_; // saves the prediction of the node
    
  bool hasChildren_;
  Node* leftChild_;
  Node* rightChild_;

};

#endif
