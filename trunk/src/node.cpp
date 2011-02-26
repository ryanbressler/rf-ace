#include<iostream>
#include<cassert>
#include "node.hpp"

Node::Node(int nsamples, bool isregr):
  isregr_(isregr),
  trainsampleics_(nsamples),
  testsampleics_(nsamples),
  ntrainsamples_(0),
  ntestsamples_(0),
  haschildren_(false)
{
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

bool Node::descend(cat_t value, int& child)
{
  if(!haschildren_)
    {
      return(false);
    }

  if(classet_.find(value) != classet_.end()) 
    {
      child = leftchild_;
    } 
  else 
    {
      child = rightchild_;
    }

  return(true);
}

bool Node::descend(num_t value, int& child)
{
  if(!haschildren_)
    {
      return(false);
    }

  if(value <= threshold_) 
    {
      child = leftchild_;
    } 
  else 
    {
      child = rightchild_;
    }

  return(true);
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

void Node::reset_testsample_ics()
{
  ntestsamples_ = 0;
}

bool Node::is_leaf()
{
  return(haschildren_);
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
	       << ") and " << threshold_ << "<x right (node " << rightchild_ << ")" << endl;
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

