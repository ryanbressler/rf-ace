#include "distributions.hpp"

#include "utils.hpp"

using datadefs::num_t;

distributions::Random::Random():
  rand_(0,datadefs::MAX_IDX) {
  this->seed( distributions::generateSeed() );
}

distributions::Random::Random(size_t seed):
  rand_(0,datadefs::MAX_IDX) {

  this->seed(seed);

}

distributions::Random::~Random() {
  
}

void distributions::Random::seed(size_t seed) {
  eng_.seed(seed);
}

size_t distributions::Random::integer() {
  return( rand_(eng_) );
}

num_t distributions::Random::uniform() {

  return( 1.0 * rand_(eng_) / ( datadefs::MAX_IDX + 1 ) );

}

distributions::PMF::PMF(const vector<num_t>& weights) {
  
  size_t n = weights.size();

  num_t sum = 0.0;
  
  for ( size_t i = 0; i < n; ++i ) {
    assert( weights[i] >= 0.0 );
    sum += weights[i];
  } 
  
  prob_.resize(n);
  alias_.resize(n);
  
  vector<size_t> HL(n);
  vector<size_t>::iterator H(HL.begin()-1);
  vector<size_t>::iterator L(HL.end());

  for ( size_t i = 0; i < n; ++i ) {
    prob_[i] = weights[i] / sum * n;
    if ( prob_[i] < 1.0 ) {
      *(++H) = i;
    } else {
      *(--L) = i;
    }
  }

  for ( size_t k = 0; k < n-1; k++ ) {
    size_t i = HL[k];
    size_t j = *L;
    alias_[i] = j;
    prob_[j] += prob_[i] - 1.0;
    if ( prob_[j] < 1.0 ) {
      L++;
    }
    if ( L >= HL.end() ) {
      break;
    }
  }
  
}

distributions::PMF::~PMF() { }

size_t distributions::PMF::sample(distributions::Random* random) const {

  size_t i = random->integer() % prob_.size(); 

  return( random->uniform() < prob_[i] ? i : alias_[i] );

}

