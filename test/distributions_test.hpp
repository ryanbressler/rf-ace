#ifndef DISTRIBUTIONSTEST_HPP
#define DISTRIBUTIONSTEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "distributions.hpp"
#include <unordered_map>

class DistributionsTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( DistributionsTest );
  CPPUNIT_TEST( test_randint );
  CPPUNIT_TEST( test_uniform );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();
  
  void test_randint();
  void test_uniform();
  
private:

  distributions::RandInt randInt_;

};

void DistributionsTest::setUp() {

  //cout << "Testing distributions.hpp: " << flush;
  randInt_.seed(0);

}

void DistributionsTest::tearDown() {

  //cout << " DONE" << endl;

}

void DistributionsTest::test_randint() {

  // Make two identical random integer generators
  distributions::RandInt randGen1(0);
  distributions::RandInt randGen2(0);

  // Test that rand1 and rand2 stay in sync
  for ( size_t i = 0; i < 1000; ++i ) {
    
    size_t r1 = randGen1();
    size_t r2 = randGen2();
    //cout << " " << r1 << "-" << r2 << " ";
    CPPUNIT_ASSERT( randGen1() == randGen2() );
  }

  unordered_map<size_t,size_t> hist;

  size_t maxIdx = 1000;

  for ( size_t i = 0; i < maxIdx; ++i ) {
    hist[i] = 0;
  }

  for ( size_t i = 0; i < 100000; ++i ) {
    //size_t r = rand1() % maxIdx;
    ++hist[ randGen1() % maxIdx ];
  }

  size_t nZeroCounts = 0;

  for ( size_t i = 0; i < maxIdx; ++i ) {
    if ( hist[i] == 0 ) ++nZeroCounts;
  }

  // We allow there to be at most two indices that never got sampled during
  // 100k random number generation rounds
  CPPUNIT_ASSERT( nZeroCounts <= 2 );

}

void DistributionsTest::test_uniform() {

  datadefs::num_t r_min = datadefs::NUM_INF;
  datadefs::num_t r_max = 0.0;

  for ( size_t i = 0; i < 100000; ++i ) {
    datadefs::num_t r = randInt_.uniform();
    CPPUNIT_ASSERT( 0.0 <= r && r <= 1.0 );

    //cout << " " << r; 

    if ( r_min > r ) r_min = r;
    if ( r_max < r ) r_max = r;
    
  }
  
  CPPUNIT_ASSERT( r_max > r_min );

  //cout << "x ~ U(" << r_min << "," << r_max << ")" << endl; 
  
  CPPUNIT_ASSERT( fabs( 1 - r_max - r_min ) < 0.0001 );

}




// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( DistributionsTest );

#endif // DISTRIBUTIONSTEST_HPP
