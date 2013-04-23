#ifndef DISTRIBUTIONS_NEWTEST_HPP
#define DISTRIBUTIONS_NEWTEST_HPP

#include <cstdlib>
#include <vector>
#include <unordered_map>

#include "utils.hpp"
#include "newtest.hpp"
#include "distributions.hpp"

using namespace std;

void distributions_newtest_random_integer();
void distributions_newtest_random_uniform();
void distributions_newtest_PMF();

void distributions_newtest() {

  newtest("integer()", &distributions_newtest_random_integer);
  newtest("uniform()", &distributions_newtest_random_uniform);
  newtest("pmf()", &distributions_newtest_PMF);

}

void distributions_newtest_random_integer() {

  // Make two identical random integer generators
  distributions::Random randGen1(0);
  distributions::Random randGen2(0);

  bool stayInSync = true;

  // Test that rand1 and rand2 stay in sync
  for ( size_t i = 0; i < 1000; ++i ) {

    size_t r1 = randGen1.integer();
    size_t r2 = randGen2.integer();

    if ( r1 != r2 ) {
      stayInSync = false;
      break;
    }
    
  }
  
  newassert( stayInSync );

  unordered_map<size_t,size_t> hist;

  size_t maxIdx = 1000;

  for ( size_t i = 0; i < maxIdx; ++i ) {
    hist[i] = 0;
  }

  for ( size_t i = 0; i < 100000; ++i ) {
    //size_t r = rand1() % maxIdx;
    ++hist[ randGen1.integer() % maxIdx ];
  }

  size_t nZeroCounts = 0;

  for ( size_t i = 0; i < maxIdx; ++i ) {
    if ( hist[i] == 0 ) ++nZeroCounts;
  }

  // We allow there to be at most two indices that never got sampled during
  // 100k random number generation rounds
  newassert( nZeroCounts <= 2 );

}

void distributions_newtest_random_uniform() {

  num_t r_min = datadefs::NUM_INF;
  num_t r_max = 0.0;

  distributions::Random random(0);

  bool stayWithinBounds = true;

  for ( size_t i = 0; i < 100000; ++i ) {
    num_t r = random.uniform();

    if ( ! (0.0 <= r && r <= 1.0) ) {
      stayWithinBounds = false;
      break;
    }

    if ( r_min > r ) r_min = r;
    if ( r_max < r ) r_max = r;
    
  }
 
  newassert( stayWithinBounds );
  newassert( r_max > r_min );
  newassert( fabs( 1 - r_max - r_min ) < 0.0001 );

}

void distributions_newtest_PMF() {

  distributions::Random random(0);

  vector<num_t> weights = {1,2,3,5,3,1,0,1e-5};

  num_t sum = math::mean(weights) * weights.size();

  distributions::PMF pmf(weights);

  vector<num_t> PMFest(8,0.0);
  
  size_t maxIter = 1e7;
  num_t incr = 1.0/maxIter;
  
  for ( size_t i = 0; i < maxIter; ++i ) {
    PMFest[ pmf.sample(&random) ] += incr;
  }
  
  for ( size_t i = 0; i < 8; ++i ) {
    newassert( fabs( PMFest[i] - weights[i] / sum ) < 0.01 );
  }
  
}

#endif
