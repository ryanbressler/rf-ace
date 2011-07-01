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
  leftChild_ = NULL;
  rightChild_ = NULL;
}

void Node::setSplitter(size_t splitter, set<num_t> classSet, Node& leftChild, Node& rightChild)
{
  assert(!hasChildren_);
  
  isSplitterNumerical_ = false;

  splitter_ = splitter;
  classSet_ = classSet;

  leftChild_ = &leftChild;
  rightChild_ = &rightChild;
  hasChildren_ = true;
}

void Node::setSplitter(size_t splitter, num_t threshold, Node& leftChild, Node& rightChild)
{
  assert(!hasChildren_);

  isSplitterNumerical_ = true;

  splitter_ = splitter;
  threshold_ = threshold;

  leftChild_ = &leftChild;
  rightChild_ = &rightChild;
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

/*
  void Node::set_impurity(num_t impurity)
  {
  impurity_ = impurity;
  }
  
  num_t Node::get_impurity()
  {
  return(impurity_);
  }
  
  void Node::reset_impurity()
  {
  impurity_ = 0.0;
  }
  
  void Node::set_prediction(num_t value)
  {
  prediction_ = value;
  }
  
  num_t Node::get_prediction()
  {
  return(prediction_);
  }
  
  void Node::set_trainidx(size_t trainidx)
  {
  trainics_.push_back(trainidx);
  }
  
  vector<size_t>* Node::get_trainics()
  {
  return(&trainics_);
  }
  
  void Node::clear_trainics()
  {
  trainics_.clear();
  }
*/

/*
  bool Node::has_children()
  {
  return(haschildren_);
  }
*/

void Node::print()
{
  if(isSplitterNumerical_) 
    {
      cout << "***NODE (numerical splitter)***" << endl;
    } 
  else 
    {
      cout << "***NODE (categorical splitter)***" << endl;
    }
  if(hasChildren_) 
    {
      cout << "-Splitter feature is " << splitter_ << endl;
      if(isSplitterNumerical_)
	{
	  cout << "-Feature value x<=" << threshold_ << " sends left (node " << &leftChild_ 
	       << ") and " << threshold_ << "<x right (node " << &rightChild_ << ")" << endl;
	}
      else
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

