#ifndef SPLITTER_HPP
#define SPLITTER_HPP

#include <set>
#include <string>

#include "datadefs.hpp"

using namespace std;
using datadefs::num_t;

class Splitter {

public:

  //Splitter();
  Splitter(const string& splitterName, num_t splitLeftLeqValue);
  Splitter(const string& splitterName, const set<num_t>& splitLeftValues, const set<num_t>& splitRightValues);
  Splitter(const Splitter& splitter);

  ~Splitter();

  bool splitsLeft(num_t testValue);
  bool splitsRight(num_t testValue);

  inline string name() { return(splitterName_); }

  /*
    void reset();
    void reset(const string& splitterName, num_t splitLeftLeqValue);
    void reset(const string& splitterName, const set<num_t>& splitLeftValues, const set<num_t>& splitRightValues);
  */

  void print();

#ifndef TEST__
private:
#endif

  enum SplitterType {NO_SPLITTER,NUMERICAL_SPLITTER,CATEGORICAL_SPLITTER};

  SplitterType splitterType_; 
  string splitterName_;

  num_t splitLeftLeqValue_;
  set<num_t> splitLeftValues_;
  set<num_t> splitRightValues_;

};


#endif
