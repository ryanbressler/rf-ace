#ifndef MATHTEST_HPP
#define MATHTEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include <cmath>
#include "math.hpp"
#include "datadefs.hpp"

using namespace std;

class MathTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( MathTest );
  CPPUNIT_TEST( test_erf );
  CPPUNIT_TEST( test_squaredError );
  CPPUNIT_TEST( test_pearsonCorrelation );
  CPPUNIT_TEST( test_ttest );
  CPPUNIT_TEST( test_mode );
  CPPUNIT_TEST( test_gamma );
  CPPUNIT_TEST( test_incrementDecrementSquaredError );
  CPPUNIT_TEST( test_incrementDecrementSquaredFrequency );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();
  
  void test_erf();
  void test_squaredError();
  void test_pearsonCorrelation();
  void test_ttest();
  void test_mode();
  void test_gamma();
  void test_incrementDecrementSquaredError();
  void test_incrementDecrementSquaredFrequency();

};

void MathTest::setUp() {}
void MathTest::tearDown() {}

void MathTest::test_erf() {

  datadefs::num_t x;
  
  x = 0.0;
  CPPUNIT_ASSERT(math::erf(x) == 0.0);
  
  x = 1.0;
  CPPUNIT_ASSERT(fabs(math::erf(x) - 0.84292558224862812465971728670410811901092529296875)
		 < datadefs::EPS);
  
  x = datadefs::NUM_NAN;
  CPPUNIT_ASSERT(datadefs::isNAN(math::erf(x)));

}

void MathTest::test_squaredError() {

  vector<datadefs::num_t> data;
  datadefs::num_t mu = -1.0;
  datadefs::num_t se = -1.0;
  //size_t nRealValues = static_cast<size_t>(-1);
  for (int i = 0; i < 50; ++i) {
    data.push_back(static_cast<datadefs::num_t>(i));
  }
  
  mu = math::mean(data);
  se = math::squaredError(data, mu);
  CPPUNIT_ASSERT(mu == 24.5);
  CPPUNIT_ASSERT(se == 10412.5);
  CPPUNIT_ASSERT( fabs( se - math::squaredError(data) ) < datadefs::EPS );
  
  data.resize(0);

  mu = math::mean(data);
  se = math::squaredError(data, mu);
  CPPUNIT_ASSERT( datadefs::isNAN(mu) );
  CPPUNIT_ASSERT( datadefs::isNAN(se) );
  CPPUNIT_ASSERT( datadefs::isNAN(math::squaredError(data)) );
  
}

void MathTest::test_ttest() {

  vector<datadefs::num_t> x,y;

  for(size_t i = 0; i < 10; ++i) {
    x.push_back(static_cast<datadefs::num_t>(i+6.0));
    y.push_back(static_cast<datadefs::num_t>(i+4.0));
  }

  // While increasing larger and larger samples, we expect
  // p-value to decrease gradually
  CPPUNIT_ASSERT( fabs( math::ttest(x,y) - 0.08 ) < 0.01 );
  CPPUNIT_ASSERT( fabs( math::ttest(y,x) - 0.92 ) < 0.01 );
  x.push_back(16.0);
  CPPUNIT_ASSERT( fabs( math::ttest(x,y) - 0.04 ) < 0.01 );
  CPPUNIT_ASSERT( fabs( math::ttest(y,x) - 0.96 ) < 0.01 );
  y.pop_back();
  CPPUNIT_ASSERT( fabs( math::ttest(x,y) - 0.02 ) < 0.01 );
  CPPUNIT_ASSERT( fabs( math::ttest(y,x) - 0.98 ) < 0.01 );
  x.push_back(17.0);
  CPPUNIT_ASSERT( fabs( math::ttest(x,y) - 0.01 ) < 0.01 );
  x.push_back(18.0);
  CPPUNIT_ASSERT( fabs( math::ttest(x,y) - 0.007 ) < 0.01 );
  x.push_back(19.0);
  CPPUNIT_ASSERT( fabs( math::ttest(x,y) - 0.005 ) < 0.01 );
  x.push_back(20.0);
  CPPUNIT_ASSERT( fabs( math::ttest(x,y) - 0.003 ) < 0.01 );
  y.pop_back();
  CPPUNIT_ASSERT( fabs( math::ttest(x,y) - 0.002 ) < 0.01 );

  datadefs::num_t p_old = 1.0;

  // While increasing values in x, we expect p-value to increase
  for ( size_t i = 0; i < 100; ++i ) {
    for ( size_t j = 0; j < x.size(); ++j ) {
      ++x[j];
    }

    datadefs::num_t p_new = math::ttest(x,y);
    //cout << "TEST " << i << ": " << p_new << " <= " << p_old << endl;
    CPPUNIT_ASSERT( p_new <= p_old );
    CPPUNIT_ASSERT( p_new < 0.002 );
    p_old = p_new;
  }

  // Taking the t-test backwards should yield a large p-value
  CPPUNIT_ASSERT( math::ttest(y,x) > 0.999 );

  p_old = 0.0;

  // While increasing values in y, we expect p-values to increase
  for ( size_t i = 0; i < 100; ++i ) {
    for ( size_t j = 0; j < y.size(); ++j ) {
      y[j] += 2;
    }

    datadefs::num_t p_new = math::ttest(x,y);
    //cout << "TEST " << i << ": " << p_new << " <= " << p_old << endl;
    CPPUNIT_ASSERT( p_new >= p_old );
    //CPPUNIT_ASSERT( p_new < 0.002 );
    p_old = p_new;
  }

  // We should now have a large p-value in p_old
  CPPUNIT_ASSERT( p_old > 0.999 );

  x.clear();
  x.resize(20,5.0);
  y.clear();
  y.resize(20,6.0);

  // p-value with constant x smaller than y should be 1.0
  CPPUNIT_ASSERT( fabs( math::ttest(x,y) - 1.0 ) < datadefs::EPS );

  y.clear();
  y.resize(20,4.0);

  // p-value with constant x larger than y should be EPS
  CPPUNIT_ASSERT( fabs( math::ttest(x,y) ) == datadefs::EPS );

  y.clear();
  y.resize(20,5.0);

  // p-value with equal constant data should be 0.5
  CPPUNIT_ASSERT( fabs( math::ttest(x,y) - 0.5 ) < datadefs::EPS );

  p_old = 0.5;

  // Start increasing x slightly to make sure p-values behave smoothly
  for ( size_t i = 0; i < 100; ++ i) {
    x[0] += 10*datadefs::EPS;

    datadefs::num_t p_new = math::ttest(x,y);

    CPPUNIT_ASSERT( p_new <= p_old );

    p_old = p_new;
  }

  x[0] = 5.0;

  p_old = 0.5;

  // Start decreasing x slightly to make sure p-values behave smoothly
  for ( size_t i = 0; i < 100; ++ i) {
    x[0] -= 10*datadefs::EPS;

    datadefs::num_t p_new = math::ttest(x,y);

    CPPUNIT_ASSERT( p_new >= p_old );

    p_old = p_new;
  }

  x.clear();
  y.clear();

  for( size_t i = 0; i < 10; ++i ) {
    x.push_back(static_cast<datadefs::num_t>(i));
    y.push_back(static_cast<datadefs::num_t>(i));
  }

  x[9] = static_cast<datadefs::num_t>(9.1);

  //cout << math::ttest(x,y) << endl;

  CPPUNIT_ASSERT( fabs( math::ttest(x,y) - 0.497 ) < 0.01 );

}

void MathTest::test_mode() {

  vector<datadefs::num_t> x(20,0);
  CPPUNIT_ASSERT( math::mode<datadefs::num_t>(x) == 0 );

  vector<size_t> y(20);
  datadefs::range(y);
  CPPUNIT_ASSERT( math::mode<size_t>(y) == 0 );

  y[19] = 18;
  CPPUNIT_ASSERT( math::mode<size_t>(y) == 18 );

  y[18] = 17;
  CPPUNIT_ASSERT( math::mode<size_t>(y) == 17 );

}

void MathTest::test_gamma() {

  vector<datadefs::num_t> data;
  size_t numClasses = 5;

  for (int i = 0; i < 50; ++i) {
    data.push_back(static_cast<datadefs::num_t>(i));
  }

  datadefs::num_t leafGamma = math::gamma(data, numClasses);
  CPPUNIT_ASSERT(leafGamma == -0.025);

}

void MathTest::test_pearsonCorrelation() {
  vector<datadefs::num_t> x;
  vector<datadefs::num_t> y;
  vector<datadefs::num_t> y2;
  datadefs::num_t corr;

  for (int i = 0; i < 50; ++i) {
    x.push_back(static_cast<datadefs::num_t>(i));
    y.push_back(static_cast<datadefs::num_t>(i));
    y2.push_back(static_cast<datadefs::num_t>(49-i));
  }
  corr = math::pearsonCorrelation(x, y);
  CPPUNIT_ASSERT(corr == 1.0);

  corr = math::pearsonCorrelation(x, y2);
  CPPUNIT_ASSERT(corr == -1.0);

  // x and y are defined as const in our signature. Since we're not
  //  testing edge cases of non-trivial memory corruption, we ignore it here.

  x.clear();
  y.clear();
  y2.clear();

  corr = math::pearsonCorrelation(x, y);
  CPPUNIT_ASSERT( datadefs::isNAN(corr) );

  for (int i = 0; i < 50; ++i) {
    x.push_back(datadefs::NUM_NAN);
    y.push_back(datadefs::NUM_NAN);
  }

  corr = math::pearsonCorrelation(x, y);
  CPPUNIT_ASSERT( datadefs::isNAN(corr) );
}


void MathTest::test_incrementDecrementSquaredError() {
  
  datadefs::num_t se = 0.0;
  datadefs::num_t mu = 0.0;
  size_t n = 0;

  // TEST: increment three values, 1.0, 2.0, 3.0
  
  ++n;
  math::incrementSquaredError(1.0,n,mu,se);
  CPPUNIT_ASSERT( fabs( mu - 1.0) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( se - 0.0) < datadefs::EPS );
  
  ++n;
  math::incrementSquaredError(2.0,n,mu,se);
  CPPUNIT_ASSERT( fabs( mu - 1.5 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( se - 0.5 ) < datadefs::EPS ); 

  ++n;
  math::incrementSquaredError(3.0,n,mu,se);
  CPPUNIT_ASSERT( fabs( mu - 2.0 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( se - 2.0 ) < datadefs::EPS );

  // TEST: decrement the same values

  --n;
  math::decrementSquaredError(3.0,n,mu,se);
  CPPUNIT_ASSERT( fabs( mu - 1.5 ) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( se - 0.5 ) < datadefs::EPS );

  --n;
  math::decrementSquaredError(2.0,n,mu,se);
  CPPUNIT_ASSERT( fabs( mu - 1.0) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( se - 0.0) < datadefs::EPS );

  --n;
  math::decrementSquaredError(1.0,n,mu,se);
  CPPUNIT_ASSERT( fabs( mu - 0.0) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( se - 0.0) < datadefs::EPS );

  // We're back at square one
  CPPUNIT_ASSERT( n == 0 );

  // TEST: increment 10 times the same value
  while( n < 10 ) {
    
    ++n;
    math::incrementSquaredError(1.5,n,mu,se);
    CPPUNIT_ASSERT( fabs( mu - 1.5) < datadefs::EPS );
    CPPUNIT_ASSERT( fabs( se - 0.0) < datadefs::EPS );

  }

  // TEST: decrement 10 times the same value
  while( n >= 2 ) {

    --n;
    math::decrementSquaredError(1.5,n,mu,se);
    CPPUNIT_ASSERT( fabs( mu - 1.5) < datadefs::EPS );
    CPPUNIT_ASSERT( fabs( se - 0.0) < datadefs::EPS );

  }

  --n;
  math::decrementSquaredError(1.5,n,mu,se);
  CPPUNIT_ASSERT( fabs( mu - 0.0) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs( se - 0.0) < datadefs::EPS );

  
}


void MathTest::test_incrementDecrementSquaredFrequency() {

  map<num_t,size_t> freq;
  size_t sqFreq = 0;


  // TEST: increment 1,1,1,2,2
  math::incrementSquaredFrequency(1,freq,sqFreq);
  CPPUNIT_ASSERT( sqFreq == 1 );

  math::incrementSquaredFrequency(1,freq,sqFreq);
  CPPUNIT_ASSERT( sqFreq == 4 );

  math::incrementSquaredFrequency(1,freq,sqFreq);
  CPPUNIT_ASSERT( sqFreq == 9 );

  math::incrementSquaredFrequency(2,freq,sqFreq);
  CPPUNIT_ASSERT( sqFreq == 10 );

  math::incrementSquaredFrequency(2,freq,sqFreq);
  CPPUNIT_ASSERT( sqFreq == 13 );

  CPPUNIT_ASSERT( freq.size() == 2 );


  // TEST: decrement 2,2,1,1,1
  math::decrementSquaredFrequency(2,freq,sqFreq);
  CPPUNIT_ASSERT( sqFreq == 10 );

  math::decrementSquaredFrequency(2,freq,sqFreq);
  CPPUNIT_ASSERT( sqFreq == 9 );

  math::decrementSquaredFrequency(1,freq,sqFreq);
  CPPUNIT_ASSERT( sqFreq == 4 );

  math::decrementSquaredFrequency(1,freq,sqFreq);
  CPPUNIT_ASSERT( sqFreq == 1 );

  math::decrementSquaredFrequency(1,freq,sqFreq);
  CPPUNIT_ASSERT( sqFreq == 0 );

  CPPUNIT_ASSERT( freq.size() == 0 );


  // TEST: increment 1,2,3,...,10
  for( size_t i = 1; i <= 10; ++i ) {
    math::incrementSquaredFrequency(static_cast<num_t>(i),freq,sqFreq);
    CPPUNIT_ASSERT( freq.size() == i );
    CPPUNIT_ASSERT( sqFreq == i );
  }


  // TEST: decrement 1,2,3,...,10
  for( size_t i = 1; i <= 10; ++i ) {
    math::decrementSquaredFrequency(static_cast<num_t>(i),freq,sqFreq);
    CPPUNIT_ASSERT( freq.size() == 10 - i );
    CPPUNIT_ASSERT( sqFreq == 10 - i );
  }


}



// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( MathTest );

#endif // MATHTEST_HPP
