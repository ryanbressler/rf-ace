#ifndef PARTITIONSEQUENCE_HPP
#define PARTITIONSEQUENCE_HPP

#include<cstdlib>
#include<vector>

using namespace std;

class PartitionSequence {

public:

  PartitionSequence();
  ~PartitionSequence();

  bool isAdded(const size_t pos);
  size_t at(const size_t pos);

  size_t upperLimit(const size_t nBits);

private:

  size_t grayCode(const size_t pos);

  

};

#endif
