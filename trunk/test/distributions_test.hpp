#ifndef DISTRIBUTIONSTEST_HPP
#define DISTRIBUTIONSTEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "distributions.hpp"
#include "math.hpp"
#include <unordered_map>

class DistributionsTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( DistributionsTest );
  CPPUNIT_TEST( test_randint );
  CPPUNIT_TEST( test_uniform );
  CPPUNIT_TEST( test_invcdf1 );
  CPPUNIT_TEST( test_invcdf2 );
  CPPUNIT_TEST( test_invcdf3 );
  CPPUNIT_TEST( test_invcdf4 );
  CPPUNIT_TEST( test_invcdf5 );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();
  
  void test_randint();
  void test_uniform();
  void test_invcdf1();
  void test_invcdf2();
  void test_invcdf3();
  void test_invcdf4();
  void test_invcdf5();

private:

  distributions::Random random_;

};

void DistributionsTest::setUp() {

  //cout << "Testing distributions.hpp: " << flush;
  random_.seed(0);

}

void DistributionsTest::tearDown() {

  //cout << " DONE" << endl;

}

void DistributionsTest::test_randint() {

  // Make two identical random integer generators
  distributions::Random randGen1(0);
  distributions::Random randGen2(0);

  // Test that rand1 and rand2 stay in sync
  for ( size_t i = 0; i < 1000; ++i ) {
    
    size_t r1 = randGen1.integer();
    size_t r2 = randGen2.integer();
    //cout << " " << r1 << "-" << r2 << " ";
    CPPUNIT_ASSERT( randGen1.integer() == randGen2.integer() );
  }

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
  CPPUNIT_ASSERT( nZeroCounts <= 2 );

}

void DistributionsTest::test_uniform() {

  datadefs::num_t r_min = datadefs::NUM_INF;
  datadefs::num_t r_max = 0.0;

  for ( size_t i = 0; i < 100000; ++i ) {
    datadefs::num_t r = random_.uniform();
    CPPUNIT_ASSERT( 0.0 <= r && r < 1.0 );

    //cout << " " << r; 

    if ( r_min > r ) r_min = r;
    if ( r_max < r ) r_max = r;
    
  }
  
  CPPUNIT_ASSERT( r_max > r_min );

  //cout << "x ~ U(" << r_min << "," << r_max << ")" << endl; 
  
  CPPUNIT_ASSERT( fabs( 1 - r_max - r_min ) < 0.0001 );

}

void DistributionsTest::test_invcdf1() {

  distributions::PMF pmf({0.2,0.2,0.2,0.2,0.2});

  CPPUNIT_ASSERT( pmf.icdf(0.00) == 0 );
  CPPUNIT_ASSERT( pmf.icdf(0.05) == 0 );
  CPPUNIT_ASSERT( pmf.icdf(0.15) == 0 );
  CPPUNIT_ASSERT( pmf.icdf(0.25) == 1 );
  CPPUNIT_ASSERT( pmf.icdf(0.35) == 1 );
  CPPUNIT_ASSERT( pmf.icdf(0.45) == 2 );
  CPPUNIT_ASSERT( pmf.icdf(0.55) == 2 );
  CPPUNIT_ASSERT( pmf.icdf(0.65) == 3 );
  CPPUNIT_ASSERT( pmf.icdf(0.75) == 3 );
  CPPUNIT_ASSERT( pmf.icdf(0.85) == 4 );
  CPPUNIT_ASSERT( pmf.icdf(0.95) == 4 );
  CPPUNIT_ASSERT( pmf.icdf(0.999999999999) == 4 );

}

void DistributionsTest::test_invcdf2() {
  
  distributions::PMF pmf({1,2,3,0,4});

  CPPUNIT_ASSERT( pmf.icdf(0.00) == 0 );
  CPPUNIT_ASSERT( pmf.icdf(0.05) == 0 );
  CPPUNIT_ASSERT( pmf.icdf(0.15) == 1 );
  CPPUNIT_ASSERT( pmf.icdf(0.25) == 1 );
  CPPUNIT_ASSERT( pmf.icdf(0.35) == 2 );
  CPPUNIT_ASSERT( pmf.icdf(0.45) == 2 );
  CPPUNIT_ASSERT( pmf.icdf(0.55) == 2 );
  CPPUNIT_ASSERT( pmf.icdf(0.65) == 4 );
  CPPUNIT_ASSERT( pmf.icdf(0.75) == 4 );
  CPPUNIT_ASSERT( pmf.icdf(0.85) == 4 );
  CPPUNIT_ASSERT( pmf.icdf(0.95) == 4 );
  CPPUNIT_ASSERT( pmf.icdf(0.99999999999) == 4 );

} 

void DistributionsTest::test_invcdf3() {

  distributions::PMF pmf({0,0,1,0,0});

  for ( size_t i = 0; i < 100; ++i ) {
    CPPUNIT_ASSERT( pmf.icdf(static_cast<datadefs::num_t>(0.01*i)) == 2 );
  }

}

void DistributionsTest::test_invcdf4() {

  vector<datadefs::num_t> weights = {1,2,3,5,3,1,0,1e-5};

  datadefs::num_t sum = math::mean(weights) * weights.size();

  distributions::PMF pmf(weights);

  vector<datadefs::num_t> PMFest(8,0.0);

  distributions::Random random;

  for ( size_t i = 0; i < 1e7; ++i ) {
    PMFest[ pmf.icdf(random.uniform()) ] += 1e-7;
  }

  for ( size_t i = 0; i < 8; ++i ) {
    CPPUNIT_ASSERT( fabs( PMFest[i] - weights[i] / sum ) < 0.01 );
  }

}

void DistributionsTest::test_invcdf5() {

  vector<datadefs::num_t> weights = {0,1};

  distributions::PMF pmf(weights);

  distributions::Random random;

  for ( size_t i = 0; i < 100; ++i ) {
    CPPUNIT_ASSERT( pmf.icdf(random.uniform()) == 1 );
  }

}

// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( DistributionsTest );

#endif // DISTRIBUTIONSTEST_HPP
