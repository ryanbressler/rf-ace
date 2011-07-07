#include<iostream>
#include<cassert>
#include "node.hpp"

Node::Node():
  isSplitterNumerical_(true),
  hasChildren_(false),
  leftChild_(NULL),
  rightChild_(NULL)
{
}

Node::~Node()
{
  if(hasChildren_)
    {
      delete leftChild_;
      delete rightChild_;
    }
}

void Node::setSplitter(size_t splitter, set<num_t> classSet)
{
  assert(!hasChildren_);
  
  isSplitterNumerical_ = false;

  splitter_ = splitter;
  classSet_ = classSet;

  leftChild_ = new Node;
  rightChild_ = new Node;
  hasChildren_ = true;
}

void Node::setSplitter(size_t splitter, num_t threshold)
{
  assert(!hasChildren_);
  isSplitterNumerical_ = true;

  splitter_ = splitter;
  threshold_ = threshold;

  leftChild_ = new Node;
  rightChild_ = new Node;
  hasChildren_ = true;
}

/*
  int Node::getSplitter()
  {
  assert(hasChildren_);
  return(splitter_);
  }
*/

Node* Node::percolateData(num_t value)
{

  if(!hasChildren_) { return(this); }

  if(isSplitterNumerical_)
    {
      if(value <= threshold_) { return(leftChild_); } else { return(rightChild_); }
    }
  else
    {
      if(classSet_.find(value) != classSet_.end()) { return(leftChild_); } else { return(rightChild_); }
    }
}

Node* Node::leftChild()
{
  if(hasChildren_)
    {
      return(leftChild_);
    }
  else
    {
      return(NULL);
    }
}

Node* Node::rightChild()
{
  if(hasChildren_)
    {
      return(rightChild_);
    }
  else
    {
      return(NULL);
    }
}

size_t Node::nNodes()
{
  size_t n = 1;
  this->recursiveNDescendantNodes(n);
  return(n);
}

void Node::recursiveNDescendantNodes(size_t& n)
{
  if(!hasChildren_)
    {
      return;
    }
  else
    {
      n += 2;
      leftChild_->recursiveNDescendantNodes(n);
      rightChild_->recursiveNDescendantNodes(n);
    }

}

void Node::setPrediction(num_t value)
{
  prediction_ = value;
}

num_t Node::getPrediction()
{
  return(prediction_);
}

/*  
  void Node::set_trainidx(size_t trainidx)
  {
  cout << "-Feature value x in {" << *classSet_.begin();
  for(set<num_t>::const_iterator it = ++classSet_.begin(); it != classSet_.end(); ++it)
  {
  cout << "," << *it;
  }
  cout << "} sends left (node " << &leftChild_ << ") and otherwise right (node " << &rightChild_ << ")" << endl;
  }
  
  } 
  else 
  {
  cout << "-Leaf node" << endl;
  }
  //cout << "-" << impurity_ << " impurity" << endl << endl;
  }
  
  void Node::print_compact()
  {
  cout << "***NODE PRINT IMPLEMENTATION MISSING***" << endl;
  }
*/

