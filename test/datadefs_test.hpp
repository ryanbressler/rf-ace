#ifndef DATADEFSTEST_HPP
#define DATADEFSTEST_HPP

#include <algorithm>
#include <limits>
#include <vector>
#include <cppunit/extensions/HelperMacros.h>
#include "datadefs.hpp"
#include "errno.hpp"

// !! TODO: Break out each individual test case instead of doing these purely
// !!  by method. Alternately, add comments visually separating each one. 
class DataDefsTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( DataDefsTest );
  CPPUNIT_TEST( test_size_tIsSigned);
  CPPUNIT_TEST( test_strv2catv );
  CPPUNIT_TEST( test_strv2numv );
  CPPUNIT_TEST( test_str2num );
  CPPUNIT_TEST( test_meanVals );
  CPPUNIT_TEST( test_mean );
  CPPUNIT_TEST( test_mode );
  CPPUNIT_TEST( test_gamma );
  CPPUNIT_TEST( test_cardinality );
  CPPUNIT_TEST( test_sqerr );
  CPPUNIT_TEST( test_countRealValues );
  CPPUNIT_TEST( test_count_freq );
  CPPUNIT_TEST( test_map_data );
  CPPUNIT_TEST( test_gini );
  CPPUNIT_TEST( test_gini );
  CPPUNIT_TEST( test_sqfreq );
  CPPUNIT_TEST( test_forward_sqfreq );
  CPPUNIT_TEST( test_forward_backward_sqfreq );
  CPPUNIT_TEST( test_range );
  CPPUNIT_TEST( test_sortDataAndMakeRef );
  CPPUNIT_TEST( test_ttest );
  CPPUNIT_TEST( test_utest );
  CPPUNIT_TEST( test_erf );
  CPPUNIT_TEST( test_regularized_betainc );
  CPPUNIT_TEST( test_spearman_correlation );
  CPPUNIT_TEST( test_pearson_correlation );
  CPPUNIT_TEST( test_percentile );
  CPPUNIT_TEST( test_isNAN );
  CPPUNIT_TEST( test_containsNAN );
  CPPUNIT_TEST( test_forward_sqerr );
  CPPUNIT_TEST( test_forward_backward_sqerr );
  CPPUNIT_TEST( test_increasingOrderOperator );
  CPPUNIT_TEST( test_decreasingOrderOperator );
  CPPUNIT_TEST( test_freqIncreasingOrderOperator );
  CPPUNIT_TEST( test_make_pairedv );
  CPPUNIT_TEST( test_separate_pairedv );
  CPPUNIT_TEST( test_sortFromRef );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();
  
  void test_size_tIsSigned();
  void test_strv2catv();
  void test_strv2numv();
  void test_str2num();
  void test_meanVals();
  void test_mean();
  void test_mode();
  void test_gamma();
  void test_cardinality();
  void test_sqerr();
  void test_countRealValues();
  void test_count_freq();
  void test_map_data();
  void test_gini();
  void test_sqfreq();
  void test_forward_sqfreq();
  void test_forward_backward_sqfreq();
  void test_range();
  void test_sortDataAndMakeRef();
  void test_ttest();
  void test_utest();
  void test_erf();
  void test_regularized_betainc();
  void test_spearman_correlation();
  void test_pearson_correlation();
  void test_percentile();
  void test_isNAN();
  void test_containsNAN();
  void test_forward_sqerr();
  void test_forward_backward_sqerr();
  void test_increasingOrderOperator();
  void test_decreasingOrderOperator();
  void test_freqIncreasingOrderOperator();
  void test_make_pairedv();
  void test_separate_pairedv();
  void test_sortFromRef();
};

void DataDefsTest::setUp() {}
void DataDefsTest::tearDown() {}

void DataDefsTest::test_size_tIsSigned() {
  
  // Ensure size_t is conformant with the 1999 ISO C standard (C99)
  CPPUNIT_ASSERT(!numeric_limits<size_t>::is_signed);
}

void DataDefsTest::test_strv2catv() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_strv2numv() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_str2num() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_meanVals() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_mean() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_mode() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_gamma() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_cardinality() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_sqerr() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_countRealValues() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_count_freq() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_map_data() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_gini() {
  // !! Note: two signatures!
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_sqfreq() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_forward_sqfreq() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_forward_backward_sqfreq() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_range() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_sortDataAndMakeRef() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_ttest() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_utest() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_erf() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_regularized_betainc() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_spearman_correlation() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_pearson_correlation() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_percentile() {
  // CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_isNAN() {

  // Test for isNAN(&string)
  CPPUNIT_ASSERT(datadefs::isNAN("NA"));
  CPPUNIT_ASSERT(datadefs::isNAN("NAN"));
  CPPUNIT_ASSERT(datadefs::isNAN("?"));

  CPPUNIT_ASSERT(!datadefs::isNAN("2"));
  CPPUNIT_ASSERT(!datadefs::isNAN("@data"));
  CPPUNIT_ASSERT(!datadefs::isNAN("NAte"));

  // Test for isNAN(num_t)
  
  CPPUNIT_ASSERT(numeric_limits<datadefs::num_t>::has_quiet_NaN);
  CPPUNIT_ASSERT(datadefs::isNAN(
                   numeric_limits<datadefs::num_t>::quiet_NaN()));
  
  if (numeric_limits<datadefs::num_t>::has_infinity) {
    CPPUNIT_ASSERT(!datadefs::isNAN(
                     numeric_limits<datadefs::num_t>::infinity()));
    CPPUNIT_ASSERT(!datadefs::isNAN(
                     -numeric_limits<datadefs::num_t>::infinity()));
  }
  
  CPPUNIT_ASSERT(!datadefs::isNAN((datadefs::num_t)0.0));
  CPPUNIT_ASSERT(!datadefs::isNAN((datadefs::num_t)-1.0));
  CPPUNIT_ASSERT(!datadefs::isNAN((datadefs::num_t)1.0));
}

void DataDefsTest::test_containsNAN() {

  vector<datadefs::num_t> justANaN;
  justANaN.push_back(numeric_limits<datadefs::num_t>::quiet_NaN());
  CPPUNIT_ASSERT(datadefs::containsNAN(justANaN));

  vector<datadefs::num_t> nanAtBack;
  nanAtBack.push_back(0.0);
  nanAtBack.push_back(1.0);
  nanAtBack.push_back(numeric_limits<datadefs::num_t>::quiet_NaN());
  CPPUNIT_ASSERT(datadefs::containsNAN(nanAtBack));

  vector<datadefs::num_t> nanAtFront;
  nanAtFront.push_back(numeric_limits<datadefs::num_t>::quiet_NaN());
  nanAtFront.push_back(-1.0);
  nanAtFront.push_back(0.0);
  CPPUNIT_ASSERT(datadefs::containsNAN(nanAtFront));

  vector<datadefs::num_t> noNaN;
  if (numeric_limits<datadefs::num_t>::has_infinity) {
    noNaN.push_back(-numeric_limits<datadefs::num_t>::infinity());
  }
  noNaN.push_back(0.0);
  noNaN.push_back(1.0);
  if (numeric_limits<datadefs::num_t>::has_infinity) {
    noNaN.push_back(numeric_limits<datadefs::num_t>::infinity());
  }
  CPPUNIT_ASSERT(!datadefs::containsNAN(noNaN));

  vector<datadefs::num_t> empty;
  CPPUNIT_ASSERT(!datadefs::containsNAN(empty));
  
}

void DataDefsTest::test_forward_sqerr() {
  datadefs::num_t x_n;
  size_t n;
  datadefs::num_t mu;
  datadefs::num_t se;

  // Ensure NaN repudiation
  //  !! TODO Eventually, this should be replaced with error code validation. 
  x_n = numeric_limits<datadefs::num_t>::quiet_NaN(); n = 0; mu = 0.0; se = 0.0;
  datadefs::forward_sqerr(x_n, n, mu, se);
  CPPUNIT_ASSERT(x_n != x_n);

  // Nothing else should have changed
  CPPUNIT_ASSERT(n == 0);
  CPPUNIT_ASSERT(mu == 0.0); 
  CPPUNIT_ASSERT(se == 0.0);

  // Ensure standard error is zeroed in the case of an unusable value for n
  x_n = 3.1415; n = 0; mu = 0.1337; se = 0.7331;
  datadefs::forward_sqerr(x_n, n, mu, se);
  CPPUNIT_ASSERT(x_n == 3.1415);
  CPPUNIT_ASSERT(n == 1); // Incremented
  CPPUNIT_ASSERT(fabs(mu - 3.1415) < datadefs::EPS); // mu += (x_n - mu) / n;
  CPPUNIT_ASSERT(se == 0.0);    // if no asserts thrown and !(n > 1), se = 0.0

  try {
    x_n = 3.1415; n = static_cast<size_t>(-1); mu = 0.1337; se = 0.7331;
    datadefs::forward_sqerr(x_n, n, mu, se);
    CPPUNIT_FAIL("datadefs::forward_sqerr didn't throw any exception; expected 'ERRNO_NUMERIC_OVERFLOW'");
  } catch (int e) {  // Ensure standard error is zeroed in the case of an
                     // unusable value for n
    CPPUNIT_ASSERT(e == ERRNO_NUMERIC_OVERFLOW);
    CPPUNIT_ASSERT(n == 0);

    // Nothing else should have changed
    CPPUNIT_ASSERT(x_n == 3.1415);
    CPPUNIT_ASSERT(mu == 0.1337);
    CPPUNIT_ASSERT(se == 0.7331);
  }

  x_n = 3.1415; n = 2; mu = 0.1337; se = 0.7331;
  datadefs::forward_sqerr(x_n, n, mu, se);
  CPPUNIT_ASSERT(x_n == 3.1415);
  CPPUNIT_ASSERT(n == 3); // Incremented
  CPPUNIT_ASSERT(fabs(mu - 1.136299999999999865707422941341064870357513427734375)
                  < datadefs::EPS); // mu += (x_n - mu) / n;
  CPPUNIT_ASSERT(fabs(se - 6.76434056000000172303998624556697905063629150390625)
                 < datadefs::EPS); // se += (x_n - mu) * (x_n - mu_old); 
}

void DataDefsTest::test_forward_backward_sqerr() {
  datadefs::num_t x_n;
  size_t n_left;
  datadefs::num_t mu_left;
  datadefs::num_t se_left;
  size_t n_right;
  datadefs::num_t mu_right;
  datadefs::num_t se_right;

  // Ensure NaN repudiation
  //  !! TODO Eventually, this should be replaced with error code validation. 
  x_n = numeric_limits<datadefs::num_t>::quiet_NaN();
   n_left = 0; mu_left = 0.0; se_left = 0.0;
   n_right = 0; mu_right = 0.0; se_right = 0.0;
   
  datadefs::forward_backward_sqerr(x_n, n_left, mu_left, se_left, n_right, mu_right, se_right);
  CPPUNIT_ASSERT(x_n != x_n);

  // Nothing else should have changed
  CPPUNIT_ASSERT(n_left == 0);
  CPPUNIT_ASSERT(mu_left == 0.0); 
  CPPUNIT_ASSERT(se_left == 0.0);
  CPPUNIT_ASSERT(n_right == 0);
  CPPUNIT_ASSERT(mu_right == 0.0); 
  CPPUNIT_ASSERT(se_right == 0.0);

  // Ensure standard error for both sides is zeroed in the case of an unusable value for n
  x_n = 3.1415;
   n_left = 0; mu_left = 0.1337; se_left = 0.7331;
   n_right = 1; mu_right = 0.9773; se_right = 0.3779;
   
  datadefs::forward_backward_sqerr(x_n, n_left, mu_left, se_left, n_right, mu_right, se_right);
  CPPUNIT_ASSERT(x_n == 3.1415);
  CPPUNIT_ASSERT(n_left == 1); // Incremented
  CPPUNIT_ASSERT(fabs(mu_left - 3.1415) < datadefs::EPS); // mu_left += (x_n - mu_left) / n_left;
  CPPUNIT_ASSERT(se_left == 0.0);  // if no asserts thrown and !(n > 1), se_left = 0.0
  CPPUNIT_ASSERT(n_right == 0); // Decremented
  CPPUNIT_ASSERT(mu_right == 0.0); // if no asserts thrown and n == 0,
                                   //  mu_right = 0.0
  CPPUNIT_ASSERT(se_right == 0.0); // if no asserts thrown and !(n > 1),
                                   //  se_right = 0.0

  x_n = 3.1415;
   n_left = 0; mu_left = 0.1337; se_left = 0.7331;
   n_right = 2; mu_right = 0.9773; se_right = 0.3779;

  datadefs::forward_backward_sqerr(x_n, n_left, mu_left, se_left, n_right, mu_right, se_right);
  CPPUNIT_ASSERT(x_n == 3.1415);
  CPPUNIT_ASSERT(n_left == 1); // Incremented
  CPPUNIT_ASSERT(fabs(mu_left - 3.1415) < datadefs::EPS); // mu_left += (x_n - mu_left) / n_left;
  CPPUNIT_ASSERT(se_left == 0.0);  // if no asserts thrown and !(n > 1), se_left = 0.0
  CPPUNIT_ASSERT(n_right == 1); // Decremented
  CPPUNIT_ASSERT(fabs(mu_right + 1.186900000000000066080474425689317286014556884765625)
                 < datadefs::EPS); // mu_right += (x_n - mu_right) / n_right;
  CPPUNIT_ASSERT(se_right == 0.0); // if no asserts thrown and !(n > 1),
                                   //  se_right = 0.0 

  try {
   x_n = 3.1415;
    n_left = static_cast<size_t>(-1); mu_left = 0.1337; se_left = 0.7331;
    n_right = 0; mu_right = 0.9773; se_right = 0.3779; // Note: n_right would
                                                        // underflow if we
                                                        // didn't bomb out on
                                                        // the overflow first.
    datadefs::forward_backward_sqerr(x_n, n_left, mu_left, se_left, n_right, mu_right, se_right);
    CPPUNIT_FAIL("datadefs::forward_backward_sqerr didn't throw any exception; expected 'ERRNO_NUMERIC_OVERFLOW'");
  } catch (int e) {  // Ensure standard error is zeroed in the case of an
                     // unusable value for n
    CPPUNIT_ASSERT(e == ERRNO_NUMERIC_OVERFLOW);
    CPPUNIT_ASSERT(n_left == 0);

    // Nothing else should have changed
    CPPUNIT_ASSERT(x_n == 3.1415);
    CPPUNIT_ASSERT(mu_left == 0.1337);
    CPPUNIT_ASSERT(se_left == 0.7331);
    CPPUNIT_ASSERT(n_right == 0);
    CPPUNIT_ASSERT(mu_right == 0.9773);
    CPPUNIT_ASSERT(se_right == 0.3779);
  }

  try {
   x_n = 3.1415;
    n_left = 1; mu_left = 0.1337; se_left = 0.7331;
    n_right = 0; mu_right = 0.9773; se_right = 0.3779; // Note: n_right would
                                                        // underflow if we
                                                        // didn't bomb out on
                                                        // the overflow first.
    datadefs::forward_backward_sqerr(x_n, n_left, mu_left, se_left, n_right, mu_right, se_right);
    CPPUNIT_FAIL("datadefs::forward_backward_sqerr didn't throw any exception; expected 'ERRNO_NUMERIC_UNDERFLOW'");
  } catch (int e) {  // Ensure standard error is zeroed in the case of an
                     // unusable value for n
    CPPUNIT_ASSERT(e == ERRNO_NUMERIC_UNDERFLOW);
    CPPUNIT_ASSERT(n_right == static_cast<size_t>(-1));

    CPPUNIT_ASSERT(x_n == 3.1415);
    CPPUNIT_ASSERT(n_left == 2); // Incremented
    CPPUNIT_ASSERT(mu_left == 0.1337);
    CPPUNIT_ASSERT(se_left == 0.7331);
    CPPUNIT_ASSERT(mu_right == 0.9773);
    CPPUNIT_ASSERT(se_right == 0.3779);
  }

  x_n = 3.1415;
   n_left = 2; mu_left = 0.1337; se_left = 0.7331;
   n_right = 4; mu_right = 0.9773; se_right = 0.3779; 

  datadefs::forward_backward_sqerr(x_n, n_left, mu_left, se_left, n_right, mu_right, se_right);
  CPPUNIT_ASSERT(x_n == 3.1415);
  CPPUNIT_ASSERT(n_left == 3); // Incremented
  CPPUNIT_ASSERT(fabs(mu_left - 1.136299999999999865707422941341064870357513427734375)
                  < datadefs::EPS); // mu_left += (x_n - mu_left) / n_left;
  CPPUNIT_ASSERT(fabs(se_left - 6.76434056000000172303998624556697905063629150390625)
                 < datadefs::EPS);  // se_left += (x_n - mu_left) * (x_n - mu_old);
  CPPUNIT_ASSERT(n_right == 3); // Decremented
  CPPUNIT_ASSERT(fabs(mu_right - 0.25589999999999990532018045996665023267269134521484375) 
                  < datadefs::EPS); // mu_right += (x_n - mu_right) / n_right;
  CPPUNIT_ASSERT(fabs(se_right + 5.8671155200000004725779945147223770618438720703125)
                 < datadefs::EPS);  // se_right += (x_n - mu_right) * (x_n - mu_old);
}

void DataDefsTest::test_increasingOrderOperator() {
  pair<datadefs::num_t, int> a(1.0, 2);
  pair<datadefs::num_t, int> b(2.0, 1);

  datadefs::increasingOrder<int> incOrder;
  CPPUNIT_ASSERT(incOrder.operator()(a,b));
  CPPUNIT_ASSERT(!incOrder.operator()(b,a));

  vector<pair<datadefs::num_t, int> > pairedTestVector;
  pairedTestVector.push_back(pair<datadefs::num_t, int>(5.0,6));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(4.0,5));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(3.0,4));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(2.0,3));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(1.0,2));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(0.0,1));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(-1.0,10));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(
                               numeric_limits<datadefs::num_t>::quiet_NaN(),0));

  sort(pairedTestVector.begin(), pairedTestVector.end(), incOrder);
  CPPUNIT_ASSERT(pairedTestVector[0].first == -1.0);
  CPPUNIT_ASSERT(pairedTestVector[0].second == 10);
  for (int i = 1; i < 7; ++i) {
    CPPUNIT_ASSERT(pairedTestVector[i].first == static_cast<double>(i-1));
    CPPUNIT_ASSERT(pairedTestVector[i].second == i);
  }
  CPPUNIT_ASSERT(pairedTestVector[7].first != pairedTestVector[7].first);
  CPPUNIT_ASSERT(pairedTestVector[7].second == 0);
}

void DataDefsTest::test_decreasingOrderOperator() {
  pair<datadefs::num_t, int> a(1.0, 2);
  pair<datadefs::num_t, int> b(2.0, 1);

  datadefs::decreasingOrder<int> decOrder;
  CPPUNIT_ASSERT(decOrder.operator()(b,a));
  CPPUNIT_ASSERT(!decOrder.operator()(a,b));

  vector<pair<datadefs::num_t, int> > pairedTestVector;
  pairedTestVector.push_back(pair<datadefs::num_t, int>(5.0,6));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(4.0,5));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(3.0,4));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(2.0,3));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(1.0,2));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(0.0,1));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(-1.0,10));
  pairedTestVector.push_back(pair<datadefs::num_t, int>(
                               numeric_limits<datadefs::num_t>::quiet_NaN(),0));

  sort(pairedTestVector.begin(), pairedTestVector.end(), decOrder);
  for (int i = 0; i > 6; --i) {
    CPPUNIT_ASSERT(pairedTestVector[i].first == static_cast<double>(5-i));
    CPPUNIT_ASSERT(pairedTestVector[i].second == 5-i);
  }
  CPPUNIT_ASSERT(pairedTestVector[6].first == -1.0);
  CPPUNIT_ASSERT(pairedTestVector[6].second == 10);
  CPPUNIT_ASSERT(pairedTestVector[7].first != pairedTestVector[7].first);
  CPPUNIT_ASSERT(pairedTestVector[7].second == 0);
}

void DataDefsTest::test_freqIncreasingOrderOperator() {
  pair<datadefs::num_t, size_t> a(1.0, 2);
  pair<datadefs::num_t, size_t> b(2.0, 1);

  datadefs::freqIncreasingOrder freqIncOrder;
  CPPUNIT_ASSERT(freqIncOrder.operator()(b,a));
  CPPUNIT_ASSERT(!freqIncOrder.operator()(a,b));

  vector<pair<datadefs::num_t, size_t> > pairedTestVector;
  pairedTestVector.push_back(pair<datadefs::num_t, size_t>(5.0,6));
  pairedTestVector.push_back(pair<datadefs::num_t, size_t>(4.0,5));
  pairedTestVector.push_back(pair<datadefs::num_t, size_t>(3.0,4));
  pairedTestVector.push_back(pair<datadefs::num_t, size_t>(2.0,3));
  pairedTestVector.push_back(pair<datadefs::num_t, size_t>(1.0,2));
  pairedTestVector.push_back(pair<datadefs::num_t, size_t>(0.0,1));
  pairedTestVector.push_back(pair<datadefs::num_t, size_t>(-1.0,10));
  pairedTestVector.push_back(pair<datadefs::num_t, size_t>(
                               numeric_limits<datadefs::num_t>::quiet_NaN(),0));

  sort(pairedTestVector.begin(), pairedTestVector.end(), freqIncOrder);
  CPPUNIT_ASSERT(pairedTestVector[0].first != pairedTestVector[0].first);
  CPPUNIT_ASSERT(pairedTestVector[0].second == 0);
  for (int i = 1; i < 7; ++i) {
    CPPUNIT_ASSERT(pairedTestVector[i].first == static_cast<double>(i-1));
    CPPUNIT_ASSERT(pairedTestVector[i].second == i);
  }
  CPPUNIT_ASSERT(pairedTestVector[7].first == -1.0);
  CPPUNIT_ASSERT(pairedTestVector[7].second == 10);
}

void DataDefsTest::test_make_pairedv() {
  vector<int> v1(50,1);
  vector<string> v2(50,"a");
  vector<pair<int,string> > p(50,pair<int,string>(0,""));

  datadefs::make_pairedv(v1,v2,p);
  for (int i = 0; i < 50; ++i) {
    CPPUNIT_ASSERT(v1[i] == 1);
    CPPUNIT_ASSERT(strcmp(v2[i].c_str(),"a") == 0);
    CPPUNIT_ASSERT(p[i].first == 1);
    CPPUNIT_ASSERT(strcmp(p[i].second.c_str(),"a") == 0);
  }
}

void DataDefsTest::test_separate_pairedv() {
  vector<pair<int,string> > p(50,pair<int,string>(1,"a"));
  vector<int> v1(50,0);
  vector<string> v2(50,"");

  datadefs::separate_pairedv(p,v1,v2);
  for (int i = 0; i < 50; ++i) {
    CPPUNIT_ASSERT(v1[i] == 1);
    CPPUNIT_ASSERT(strcmp(v2[i].c_str(),"a") == 0);
    CPPUNIT_ASSERT(p[i].first == 1);
    CPPUNIT_ASSERT(strcmp(p[i].second.c_str(),"a") == 0);
  }
}

void DataDefsTest::test_sortFromRef() {
  vector<int> data(50,0);
  vector<size_t> refIcs(50,0);
  for (int i = 0; i < 50; ++i) {
    data[i] = i;
    refIcs[i] = 49-i;
  }
  
  datadefs::sortFromRef<int>(data,refIcs);
  for (int i = 0; i < 50; ++i) {
    CPPUNIT_ASSERT(data[i] == 49-i);
  }
}

// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( DataDefsTest );

#endif // DATADEFSTEST_HPP
