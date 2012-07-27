#include "distributions.hpp"

using datadefs::num_t;

distributions::RandInt::RandInt():
  rand_(0,datadefs::MAX_IDX) {
  this->seed( distributions::generateSeed() );
}

distributions::RandInt::RandInt(size_t seed):
  rand_(0,datadefs::MAX_IDX) {

  this->seed(seed);

}

distributions::RandInt::~RandInt() {
  
}

void distributions::RandInt::seed(size_t seed) {
  eng_.seed(seed);
}

num_t distributions::RandInt::uniform() {

  return( 1.0 * rand_(eng_) / datadefs::MAX_IDX );

}

