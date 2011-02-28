#include<iostream>
#include<cassert>
#include "node.hpp"

Node::Node(int nsamples):
  isnum_(true),
  impurity_(0.0),
  haschildren_(false),
  leftchild_(NULL),
  rightchild_(NULL)
{
}

Node::~Node()
{
  leftchild_ = NULL;
  rightchild_ = NULL;
}

void Node::set_splitter(int splitter, set<cat_t> classet, Node& leftchild, Node& rightchild)
{
  assert(!haschildren_);
  
  isnum_ = false;

  splitter_ = splitter;
  classet_ = classet;

  leftchild_ = &leftchild;
  rightchild_ = &rightchild;
  haschildren_ = true;
}

void Node::set_splitter(int splitter, num_t threshold, Node& leftchild, Node& rightchild)
{
  assert(!haschildren_);

  isnum_ = true;

  splitter_ = splitter;
  threshold_ = threshold;

  leftchild_ = &leftchild;
  rightchild_ = &rightchild;
  haschildren_ = true;
}

int Node::get_splitter()
{
  assert(haschildren_);

  return(splitter_);
}

Node* Node::percolate(cat_t value)
{

  if(!haschildren_)
    {
      return(this);
    }

  if(classet_.find(value) != classet_.end()) 
    {
      return(leftchild_);
    } 
  else 
    {
      return(rightchild_);
    }
}

Node* Node::percolate(num_t value)
{
  
  if(!haschildren_)
    {
      return(this);
    }

  if(value <= threshold_) 
    {
      return(leftchild_);
    } 
  else 
    {
      return(rightchild_);
    }
}

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

bool Node::has_children()
{
  return(haschildren_);
}

void Node::print()
{
  if(isnum_) 
    {
      cout << "***NODE (numerical splitter)***" << endl;
    } 
  else 
    {
      cout << "***NODE (categorical splitter)***" << endl;
    }
  if(haschildren_) 
    {
      cout << "-Splitter feature is " << splitter_ << endl;
      if(isnum_)
	{
	  cout << "-Feature value x<=" << threshold_ << " sends left (node " << &leftchild_ 
	       << ") and " << threshold_ << "<x right (node " << &rightchild_ << ")" << endl;
	}
      else
	{
	  cout << "-Feature value x in {" << *classet_.begin();
	  for(set<cat_t>::const_iterator it = ++classet_.begin(); it != classet_.end(); ++it)
	    {
	      cout << "," << *it;
	    }
	  cout << "} sends left (node " << &leftchild_ << ") and otherwise right (node " << &rightchild_ << ")" << endl;
	}

    } 
  else 
    {
      cout << "-Leaf node" << endl;
    }
  cout << "-" << impurity_ << " impurity" << endl << endl;
}

void Node::print_compact()
{
  cout << "***NODE PRINT IMPLEMENTATION MISSING***" << endl;
}

