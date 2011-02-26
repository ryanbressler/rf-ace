#include<iostream>
#include<cassert>
#include "node.hpp"

Node::Node(int nsamples, bool isregr):
  isregr_(isregr),
  ntrainsamples_(0),
  ntestsamples_(0),
  haschildren_(false)
{
  if(isregr_)
    {
      vector<num_t> num_trainsamples_(nsamples);
      vector<num_t> num_testsamples_(nsamples);
    }
  else
    {
      vector<cat_t> cat_trainsamples_(nsamples);
      vector<cat_t> cat_testsamples_(nsamples);
    }
}

Node::~Node()
{

}

void Node::set_splitter(int splitter, set<cat_t> classet, int leftchild, int rightchild)
{
  assert(!isregr_);
  assert(!haschildren_);

  splitter_ = splitter;
  classet_ = classet;

  leftchild_ = leftchild;
  rightchild_ = rightchild;
  haschildren_ = true;
}

void Node::set_splitter(int splitter, num_t threshold, int leftchild, int rightchild)
{
  assert(isregr_);
  assert(!haschildren_);

  splitter_ = splitter;
  threshold_ = threshold;

  leftchild_ = leftchild;
  rightchild_ = rightchild;
  haschildren_ = true;
}

int Node::get_splitter()
{
  assert(haschildren_);

  return(splitter_);
}

int Node::percolate(cat_t value)
{
  if(classet_.find(value) != classet_.end()) {return(leftchild_);} else {return(rightchild_);}
}

int Node::percolate(num_t value)
{
  if(value <= threshold_) {return(leftchild_);} else {return(rightchild_);}
}

void Node::print()
{
  cout << endl << "***NODE PRINT***" << endl;
  if(isregr_) 
    {
      cout << "-Regression node" << endl;
    } 
  else 
    {
      cout << "-Classification node" << endl;
    }
  if(haschildren_) 
    {
      cout << "-Splitter feature is " << splitter_ << endl;
      if(isregr_)
	{
	  cout << "-Feature value x<=" << threshold_ << " sends left (node " << leftchild_ 
	       << ") and x>" << threshold_ << " right (node " << rightchild_ << ")" << endl;
	}
      else
	{
	  cout << "-Feature value x in {" << *classet_.begin();
	  for(set<cat_t>::const_iterator it = ++classet_.begin(); it != classet_.end(); ++it)
	    {
	      cout << "," << *it;
	    }
	  cout << "} sends left (node " << leftchild_ << ") and otherwise right (node " << rightchild_ << ")" << endl;
	}

    } 
  else 
    {
      cout << "-Leaf node" << endl;
    }
  
  cout << endl;
}

void Node::print_compact()
{
  cout << "***NODE PRINT IMPLEMENTATION MISSING***" << endl;
}

