#include<iostream>
#include<cassert>
#include "node.hpp"

Node::Node(int nsamples):
  isnum_(true),
  trainsampleics_(nsamples),
  testsampleics_(nsamples),
  ntrainsamples_(0),
  ntestsamples_(0),
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

Node* Node::descend(cat_t value)
{

  if(!haschildren_)
    {
      return(NULL);
    }

  if(classet_.find(value) != classet_.end()) 
    {
      return(leftchild_);
    } 
  else 
    {
      return(rightchild_);
    }

  //return(true);
}

Node* Node::descend(num_t value)
{
  
  if(!haschildren_)
    {
      return(NULL);
    }

  if(value <= threshold_) 
    {
      return(leftchild_);
    } 
  else 
    {
      return(rightchild_);
    }

  //return(true);
}

void Node::add_trainsample_idx(int idx)
{
  trainsampleics_[ntrainsamples_] = idx;
  ++ntrainsamples_;
}

void Node::add_testsample_idx(int idx)
{
  testsampleics_[ntestsamples_] = idx;
  ++ntestsamples_;
}

void Node::reset_trainsample_ics()
{
  ntrainsamples_ = 0;
}

void Node::reset_testsample_ics()
{
  ntestsamples_ = 0;
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
  cout << "-" << ntrainsamples_ << " train sample ics";
  if(ntrainsamples_)
    {
      cout << ": " << trainsampleics_[0]; 
      size_t i = 1;
      for(i = 1; i < 10 && i < ntrainsamples_; ++i)
	{
	  cout << "," << trainsampleics_[i];
	}
      if(i < ntrainsamples_)
	{
	  cout << "...";
	}
    }
  cout << endl;
  cout << "-" << ntestsamples_ << " test sample ics";
  if(ntestsamples_)
    {
      cout << ": " << testsampleics_[0];
      size_t i = 1;
      for(i = 1; i < 10 && i < ntestsamples_; ++i)
        {
          cout << "," << testsampleics_[i];
        }
      if(i < ntestsamples_)
        {
          cout << "...";
        }
    }

  cout << endl << endl;
}

void Node::print_compact()
{
  cout << "***NODE PRINT IMPLEMENTATION MISSING***" << endl;
}

