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

distributions::InvCDF::InvCDF(const vector<num_t>& weights) {

  num_t sum = 0.0;

  for ( size_t i = 0; i < weights.size(); ++i ) {
    assert( weights[i] >= 0.0 );
    sum += weights[i];
  }

  num_t cumProb = 0.0;

  assert( sum > 0.0 );

  for ( size_t i = 0; i < weights.size(); ++i ) {
    if ( weights[i] > 0.0 ) {
      cumProb += weights[i] / sum;
      icdf_[cumProb] = i;
    }
  }

}

distributions::InvCDF::~InvCDF() {
  
}

size_t distributions::InvCDF::at(const num_t prob) {
  assert(prob >= 0.0 && prob < 1.0);

  return( icdf_.upper_bound(prob)->second );
}
