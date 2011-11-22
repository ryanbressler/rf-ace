#ifndef SPLITTER_HPP
#define SPLITTER_HPP

#include <set>

#include "datadefs.hpp"

using namespace std;
using datadefs::num_t;

class Splitter {

public:

  Splitter();
  Splitter(num_t splitLeftLeqValue);
  Splitter(const set<num_t>& splitLeftValues, const set<num_t>& splitRightValues);
  Splitter(const Splitter& splitter);

  ~Splitter();

  bool splitsLeft(num_t testValue);
  bool splitsRight(num_t testValue);

  void reset();
  void reset(num_t splitLeftLeqValue);
  void reset(const set<num_t>& splitLeftValues, const set<num_t>& splitRightValues);

  void print();

#ifndef TEST__
private:
#endif

  enum SplitterType {NO_SPLITTER,NUMERICAL_SPLITTER,CATEGORICAL_SPLITTER};

  SplitterType splitterType_; 

  num_t splitLeftLeqValue_;
  set<num_t> splitLeftValues_;
  set<num_t> splitRightValues_;

};


#endif
