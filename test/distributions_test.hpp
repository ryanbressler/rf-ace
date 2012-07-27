#ifndef DISTRIBUTIONSTEST_HPP
#define DISTRIBUTIONSTEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "distributions.hpp"

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
  distributions::RandInt rand1(0);
  distributions::RandInt rand2(0);
  
  // Test that rand1 and rand2 stay in sync
  for ( size_t i = 0; i < 1000; ++i ) {
    CPPUNIT_ASSERT( rand1() == rand2() );
  }

}

void DistributionsTest::test_uniform() {

  datadefs::num_t r_min = datadefs::NUM_INF;
  datadefs::num_t r_max = 0.0;

  for ( size_t i = 0; i < 100000; ++i ) {
    datadefs::num_t r = randInt_.uniform();
    CPPUNIT_ASSERT( 0.0 <= r && r < 1.0 );

    if ( r_min > r ) r_min = r;
    if ( r_max < r ) r_max = r;
    
  }
  
  CPPUNIT_ASSERT( r_max > r_min );
  
  CPPUNIT_ASSERT( fabs( 1 - r_max - r_min ) < 0.0001 );

}


// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( DistributionsTest );

#endif // DISTRIBUTIONSTEST_HPP
