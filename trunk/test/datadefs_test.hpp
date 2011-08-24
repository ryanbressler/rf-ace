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
// !! TODO: Testing generification. This is done procedurally to be very
// !!  explicit about what is tested for now, but this can be greatly improved
// !!  with the correct function constructs.
class DataDefsTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( DataDefsTest );
  CPPUNIT_TEST( test_size_tIsSigned);
  CPPUNIT_TEST( test_strv2catv );
  CPPUNIT_TEST( test_strv2numv );
  CPPUNIT_TEST( test_str2num );
  CPPUNIT_TEST( test_mean );
  //CPPUNIT_TEST( test_mode );
  //CPPUNIT_TEST( test_gamma );
  CPPUNIT_TEST( test_cardinality );
  CPPUNIT_TEST( test_sqerr );
  CPPUNIT_TEST( test_countRealValues );
  CPPUNIT_TEST( test_count_freq );
  CPPUNIT_TEST( test_map_data );
  CPPUNIT_TEST( test_gini );
  CPPUNIT_TEST( test_sqfreq );
  CPPUNIT_TEST( test_forward_sqfreq );
  CPPUNIT_TEST( test_forward_backward_sqfreq );
  CPPUNIT_TEST( test_range );
  CPPUNIT_TEST( test_sortDataAndMakeRef );
  CPPUNIT_TEST( test_utest );
  CPPUNIT_TEST( test_erf );
  CPPUNIT_TEST( test_pearson_correlation );
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

  // Tests not yet implemented
  //CPPUNIT_TEST( test_spearman_correlation );

  // Tests for deprecated methods
  //CPPUNIT_TEST( test_ttest );
  //CPPUNIT_TEST( test_regularized_betainc );
  //CPPUNIT_TEST( test_percentile );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();
  
  void test_size_tIsSigned();
  void test_strv2catv();
  void test_strv2numv();
  void test_str2num();
  void test_mean();
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
  vector<string> strvec(51,"");
  vector<datadefs::num_t> catvec(51,0.0);
  map<string,datadefs::num_t> mapping;
  map<datadefs::num_t,string> backMapping;
  strvec[0] = "a";
  strvec[1] = "b";
  strvec[2] = "c";
  strvec[3] = "A";
  strvec[50] = "NaN";
  datadefs::strv2catv(strvec, catvec, mapping, backMapping);

  CPPUNIT_ASSERT(catvec[0] == 0.0);
  CPPUNIT_ASSERT(catvec[1] == 1.0);
  CPPUNIT_ASSERT(catvec[2] == 2.0);
  CPPUNIT_ASSERT(catvec[3] == 0.0);
  for (int i = 4; i < 50; ++i) {
    CPPUNIT_ASSERT(catvec[i] == 3.0);
  }
  CPPUNIT_ASSERT(datadefs::isNAN(catvec[50]));
}

void DataDefsTest::test_strv2numv() {
  vector<string> strvec(51,"3.0");
  vector<datadefs::num_t> catvec(51,0.0);
  strvec[0] = "0.0";
  strvec[1] = "1.0";
  strvec[2] = "2.0";
  strvec[3] = "0.00";
  strvec[50] = "NaN";
  datadefs::strv2numv(strvec, catvec);

  CPPUNIT_ASSERT(catvec[0] == 0.0);
  CPPUNIT_ASSERT(catvec[1] == 1.0);
  CPPUNIT_ASSERT(catvec[2] == 2.0);
  CPPUNIT_ASSERT(catvec[3] == 0.0);
  for (int i = 4; i < 50; ++i) {
    CPPUNIT_ASSERT(catvec[i] == 3.0);
  }
  CPPUNIT_ASSERT(datadefs::isNAN(catvec[50]));
}

void DataDefsTest::test_str2num() {
  string a("0.0");
  string b("1.0");
  string c("-1.0");
  string d("-1.0e10");

  CPPUNIT_ASSERT(datadefs::str2num(a) == 0.0);
  CPPUNIT_ASSERT(datadefs::str2num(b) == 1.0);
  CPPUNIT_ASSERT(datadefs::str2num(c) == -1.0);
  CPPUNIT_ASSERT(datadefs::str2num(d) == -1.0e10);

  // These trigger behavior that currently throws a warning, but should throw a
  //  runtime exception. !! TODO: these should become failure cases that raise,
  //  at best, exceptions for the system to handle.
  // string e("-1.0AAAA");
  // string f("");
  // string g("ABCDE");
  // string h("NaN");
  // CPPUNIT_ASSERT(datadefs::str2num(e) == -1.0);
  // CPPUNIT_ASSERT(datadefs::str2num(f) == 0.0);
  // CPPUNIT_ASSERT(datadefs::str2num(g) == 0.0);
  // CPPUNIT_ASSERT(datadefs::str2num(h) == 0.0);
}

void DataDefsTest::test_mean() {
  vector<datadefs::num_t> data;
  datadefs::num_t mu;
  size_t nRealValues;

  for (int i = 0; i < 50; ++i) {
    data.push_back(static_cast<datadefs::num_t>(i));
  }

  // Spuriously assign the values of mu and nRealValues, since they'll be
  //  flattened during function invocation
  mu = -1.0; 
  nRealValues = static_cast<size_t>(-1);
  
  datadefs::mean(data, mu, nRealValues);
  CPPUNIT_ASSERT(mu == 24.5);
  CPPUNIT_ASSERT(nRealValues == 50);

  // Our data vector is defined as const in our signature. Since we're not
  //  testing edge cases of non-trivial memory corruption, we ignore it here.

  // Interleave the original input with NaNs; verify we get the same results
  for (int i = 0; i < 50; ++i) {
    data.insert(data.begin() + (i*2), datadefs::NUM_NAN);
  }

  mu = -1.0; 
  nRealValues = static_cast<size_t>(-1);

  datadefs::mean(data, mu, nRealValues);
  CPPUNIT_ASSERT(mu == 24.5);
  CPPUNIT_ASSERT(nRealValues == 50);
  
  // Ensure a vector containing only NaNs is handled properly
  data.clear();
  for (int i = 0; i < 50; ++i) {
    data.push_back(datadefs::NUM_NAN);
  }

  mu = -1.0; 
  nRealValues = static_cast<size_t>(-1);

  datadefs::mean(data, mu, nRealValues);
  CPPUNIT_ASSERT(mu == 0.0);
  CPPUNIT_ASSERT(nRealValues == 0);
}

void DataDefsTest::test_cardinality() {
  vector<datadefs::num_t> data;
  size_t cardinality = static_cast<size_t>(-1);
  for (int i = 0; i < 50; ++i) {
    data.push_back(static_cast<datadefs::num_t>(i));
  }
  
  datadefs::cardinality(data, cardinality);
  CPPUNIT_ASSERT(cardinality == 50);

  // Interleave the original input with NaNs; verify we get the same results
  for (int i = 0; i < 50; ++i) {
    data.insert(data.begin() + (i*2), datadefs::NUM_NAN);
  }

  cardinality = static_cast<size_t>(-1);

  datadefs::cardinality(data, cardinality);
  CPPUNIT_ASSERT(cardinality == 50); 
  
  // Ensure a vector containing only NaNs is handled as expected
  data.clear();
  for (int i = 0; i < 50; ++i) {
    data.push_back(datadefs::NUM_NAN);
  }

  cardinality = static_cast<size_t>(-1); 

  datadefs::cardinality(data, cardinality);
  CPPUNIT_ASSERT(cardinality == 0);
}

void DataDefsTest::test_sqerr() {
  vector<datadefs::num_t> data;
  datadefs::num_t mu = -1.0;
  datadefs::num_t se = -1.0;
  size_t nRealValues = static_cast<size_t>(-1);
  for (int i = 0; i < 50; ++i) {
    data.push_back(static_cast<datadefs::num_t>(i));
  }
  
  datadefs::sqerr(data, mu, se, nRealValues);
  CPPUNIT_ASSERT(mu == 24.5);
  CPPUNIT_ASSERT(se == 10412.5);
  CPPUNIT_ASSERT(nRealValues == 50);
  
  // Our data vector is defined as const in our signature. Since we're not
  //  testing edge cases of non-trivial memory corruption, we ignore it here.
  
  // Interleave the original input with NaNs; verify we get the same results
  for (int i = 0; i < 50; ++i) {
    data.insert(data.begin() + (i*2), datadefs::NUM_NAN);
  }

  nRealValues = static_cast<size_t>(-1);
  
  datadefs::sqerr(data, mu, se, nRealValues);
  CPPUNIT_ASSERT(mu == 24.5);
  CPPUNIT_ASSERT(se == 10412.5);
  CPPUNIT_ASSERT(nRealValues == 50); 
  
  // Ensure a vector containing only NaNs is handled as expected
  data.clear();
  for (int i = 0; i < 50; ++i) {
    data.push_back(datadefs::NUM_NAN);
  }

  nRealValues = static_cast<size_t>(-1); 
  
  datadefs::sqerr(data, mu, se, nRealValues);
  CPPUNIT_ASSERT(mu == 0.0);
  CPPUNIT_ASSERT(se == 0.0);
  CPPUNIT_ASSERT(nRealValues == 0);
}

void DataDefsTest::test_countRealValues() {
  vector<datadefs::num_t> data;
  size_t nRealValues = static_cast<size_t>(-1);
  for (int i = 0; i < 50; ++i) {
    data.push_back(static_cast<datadefs::num_t>(i));
  }
  
  datadefs::countRealValues(data, nRealValues);
  CPPUNIT_ASSERT(nRealValues == 50);
  
  // Our data vector is defined as const in our signature. Since we're not
  //  testing edge cases of non-trivial memory corruption, we ignore it here.
  
  // Interleave the original input with NaNs; verify we get the same results
  for (int i = 0; i < 50; ++i) {
    data.insert(data.begin() + (i*2), datadefs::NUM_NAN);
  }

  nRealValues = static_cast<size_t>(-1);

  datadefs::countRealValues(data, nRealValues);
  CPPUNIT_ASSERT(nRealValues == 50); 
  
  // Ensure a vector containing only NaNs is handled as expected
  data.clear();
  for (int i = 0; i < 50; ++i) {
    data.push_back(datadefs::NUM_NAN);
  }

  nRealValues = static_cast<size_t>(-1); 
  
  datadefs::countRealValues(data, nRealValues);
  CPPUNIT_ASSERT(nRealValues == 0);
}

// !! TODO: make this test slightly more robust.
void DataDefsTest::test_count_freq() {
  vector<datadefs::num_t> data;
  map<datadefs::num_t,size_t> cat2freq;
  size_t nRealValues = static_cast<size_t>(-1);
  for (int i = 0; i < 50; ++i) {
    data.push_back(static_cast<datadefs::num_t>(i));
  }
  data.push_back(0.0);
  
  datadefs::count_freq(data, cat2freq, nRealValues);
  CPPUNIT_ASSERT(cat2freq.size() == 50);
  CPPUNIT_ASSERT(nRealValues == 51);
  
  // Our data vector is defined as const in our signature. Since we're not
  //  testing edge cases of non-trivial memory corruption, we ignore it here.
  
  // Interleave the original input with NaNs; verify we get the same results
  for (int i = 0; i < 50; ++i) {
    data.insert(data.begin() + (i*2), datadefs::NUM_NAN);
  }

  cat2freq.clear();
  nRealValues = static_cast<size_t>(-1);

  datadefs::count_freq(data, cat2freq, nRealValues);
  CPPUNIT_ASSERT(cat2freq.size() == 50);
  CPPUNIT_ASSERT(nRealValues == 51);
  
  // Ensure a vector containing only NaNs is handled as expected
  data.clear();
  for (int i = 0; i < 50; ++i) {
    data.push_back(datadefs::NUM_NAN);
  }

  cat2freq.clear();
  nRealValues = static_cast<size_t>(-1); 
  
  datadefs::count_freq(data, cat2freq, nRealValues);
  CPPUNIT_ASSERT(cat2freq.size() == 0);
  CPPUNIT_ASSERT(nRealValues == 0);
}

// !! TODO: make this test more robust
void DataDefsTest::test_map_data() {
  vector<datadefs::num_t> data;
  map<datadefs::num_t,vector<size_t> > datamap; 
  size_t nRealValues;
  for (int i = 0; i < 50; ++i) {
    data.push_back(static_cast<datadefs::num_t>(i));
  }
  data.push_back(0.0);

  datadefs::map_data(data, datamap, nRealValues);
  CPPUNIT_ASSERT(datamap.size() == 50);
  CPPUNIT_ASSERT(nRealValues == 51);
    
  // Our data vector is defined as const in our signature. Since we're not
  //  testing edge cases of non-trivial memory corruption, we ignore it here.
  
  // Interleave the original input with NaNs; verify we get the same results
  for (int i = 0; i < 50; ++i) {
    data.insert(data.begin() + (i*2), datadefs::NUM_NAN);
  }

  nRealValues = static_cast<size_t>(-1);
  datadefs::map_data(data, datamap, nRealValues);
  CPPUNIT_ASSERT(datamap.size() == 50);
  CPPUNIT_ASSERT(nRealValues == 51);
  
  // Ensure a vector containing only NaNs is handled as expected
  data.clear();
  for (int i = 0; i < 50; ++i) {
    data.push_back(datadefs::NUM_NAN);
  }

  nRealValues = static_cast<size_t>(-1); 
  datadefs::map_data(data, datamap, nRealValues);
  CPPUNIT_ASSERT(datamap.size() == 0);
  CPPUNIT_ASSERT(nRealValues == 0);
}

void DataDefsTest::test_gini() {
  vector<datadefs::num_t> data;
  datadefs::num_t giniIndex; // I dream of gini.
                             //  </obligatoryFunnyComment>
  size_t nRealValues;
  map<datadefs::num_t,size_t> cat2freq;

  for (int i = 0; i < 50; ++i) {
    data.push_back(static_cast<datadefs::num_t>(i));
  }
  data.push_back(0.0);

  datadefs::gini(data, giniIndex, nRealValues);
  CPPUNIT_ASSERT(fabs(giniIndex - 0.9796232218377547429355445274268276989459991455078125)
                  < datadefs::EPS);
  CPPUNIT_ASSERT(nRealValues == 51);
  
  datadefs::count_freq(data, cat2freq, nRealValues);
  CPPUNIT_ASSERT(cat2freq.size() == 50);
  CPPUNIT_ASSERT(nRealValues == 51);

  datadefs::gini(cat2freq, giniIndex);
  CPPUNIT_ASSERT(fabs(giniIndex - 0.9796232218377547429355445274268276989459991455078125)
                  < datadefs::EPS);
  
  // Our data vector is defined as const in our signature. Since we're not
  //  testing edge cases of non-trivial memory corruption, we ignore it here.
  
  // Interleave the original input with NaNs; verify we get the same results
  for (int i = 0; i < 50; ++i) {
    data.insert(data.begin() + (i*2), datadefs::NUM_NAN);
  }

  nRealValues = static_cast<size_t>(-1);
  datadefs::gini(data, giniIndex, nRealValues);
  CPPUNIT_ASSERT(fabs(giniIndex - 0.9796232218377547429355445274268276989459991455078125)
                  < datadefs::EPS);
  CPPUNIT_ASSERT(nRealValues == 51);
  
  datadefs::count_freq(data, cat2freq, nRealValues);
  CPPUNIT_ASSERT(cat2freq.size() == 50);
  CPPUNIT_ASSERT(nRealValues == 51);

  datadefs::gini(cat2freq, giniIndex);
  CPPUNIT_ASSERT(fabs(giniIndex - 0.9796232218377547429355445274268276989459991455078125)
                  < datadefs::EPS);
  
  // Ensure a vector containing only NaNs is handled as expected
  data.clear();
  for (int i = 0; i < 50; ++i) {
    data.push_back(datadefs::NUM_NAN);
  }

  nRealValues = static_cast<size_t>(-1); 
  datadefs::gini(data, giniIndex, nRealValues);
  CPPUNIT_ASSERT(giniIndex == 0.0);
  CPPUNIT_ASSERT(nRealValues == 0);
  
  datadefs::count_freq(data, cat2freq, nRealValues);
  CPPUNIT_ASSERT(cat2freq.size() == 0);
  CPPUNIT_ASSERT(nRealValues == 0);

  datadefs::gini(cat2freq, giniIndex);
  CPPUNIT_ASSERT(giniIndex == 0.0);

  
}

void DataDefsTest::test_sqfreq() {
  vector<datadefs::num_t> data; 
  map<datadefs::num_t,size_t> freq;
  size_t sqFreq;
  size_t nRealValues;
  for (int i = 0; i < 50; ++i) {
    data.push_back(static_cast<datadefs::num_t>(i));
  }
  data.push_back(0.0);
  
  datadefs::sqfreq(data, freq, sqFreq, nRealValues);
  CPPUNIT_ASSERT(freq.size() == 50);
  CPPUNIT_ASSERT(sqFreq == 53);
  CPPUNIT_ASSERT(nRealValues == 51);

  // Our data vector is defined as const in our signature. Since we're not
  //  testing edge cases of non-trivial memory corruption, we ignore it here.

  // Interleave the original input with NaNs; verify we get the same results
  for (int i = 0; i < 50; ++i) {
    data.insert(data.begin() + (i*2), datadefs::NUM_NAN);
  }

  freq.clear();
  sqFreq = static_cast<size_t>(-1);
  nRealValues = static_cast<size_t>(-1);

  datadefs::sqfreq(data, freq, sqFreq, nRealValues);
  CPPUNIT_ASSERT(freq.size() == 50);
  CPPUNIT_ASSERT(sqFreq == 53);
  CPPUNIT_ASSERT(nRealValues == 51);
  
  // Ensure a vector containing only NaNs is handled properly
  data.clear();
  for (int i = 0; i < 50; ++i) {
    data.push_back(datadefs::NUM_NAN);
  }
  
  freq.clear();
  sqFreq = static_cast<size_t>(-1);
  nRealValues = static_cast<size_t>(-1);
  
  datadefs::sqfreq(data, freq, sqFreq, nRealValues);
  CPPUNIT_ASSERT(freq.size() == 0);
  CPPUNIT_ASSERT(sqFreq == 0);
  CPPUNIT_ASSERT(nRealValues == 0);
}

void DataDefsTest::test_forward_sqfreq() {
  
  vector<datadefs::num_t> data; 
  map<datadefs::num_t,size_t> freq;
  size_t sqFreq;
  size_t nRealValues;
  datadefs::num_t x_n;
  size_t n;
  for (int i = 0; i < 50; ++i) {
    data.push_back(static_cast<datadefs::num_t>(i));
  }
  data.push_back(0.0);
  
  x_n = 50; n = 50;
  datadefs::sqfreq(data, freq, sqFreq, nRealValues);
  datadefs::forward_sqfreq(x_n, n, freq, sqFreq);
  CPPUNIT_ASSERT(n == 51);
  CPPUNIT_ASSERT(freq.size() == 51);
  CPPUNIT_ASSERT(sqFreq == 54);

  // Ensure NaN repudiation
  //  !! TODO Eventually, this should be replaced with error code validation. 
  x_n = datadefs::NUM_NAN; n = 0;
  datadefs::forward_sqfreq(x_n, n, freq, sqFreq);

  // Nothing should have changed
  CPPUNIT_ASSERT(n == 0);
  CPPUNIT_ASSERT(freq.size() == 51);
  CPPUNIT_ASSERT(sqFreq == 54);
  
  // Ensure standard error is zeroed in the case of an unusable value for n
  x_n = 3.1415; n = 0; 
  datadefs::forward_sqfreq(x_n, n, freq, sqFreq);
  CPPUNIT_ASSERT(n == 1); // Incremented
  CPPUNIT_ASSERT(freq.size() == 52);
  CPPUNIT_ASSERT(sqFreq == 55);
  
  // Ensure standard error is zeroed in the case of an unusable value for n
  x_n = 10; n = 0; 
  datadefs::forward_sqfreq(x_n, n, freq, sqFreq);
  CPPUNIT_ASSERT(n == 1); // Incremented
  CPPUNIT_ASSERT(freq.size() == 52);
  CPPUNIT_ASSERT(sqFreq == 58);

  try {
    x_n = 3.1415; n = static_cast<size_t>(-1); 
    datadefs::forward_sqfreq(x_n, n, freq, sqFreq);
    CPPUNIT_FAIL("datadefs::forward_sqfreq didn't throw any exception; expected 'ERRNO_NUMERIC_OVERFLOW'");
  } catch (int e) {  // Ensure standard error is zeroed in the case of an
                     // unusable value for n
    CPPUNIT_ASSERT(e == ERRNO_NUMERIC_OVERFLOW);
    CPPUNIT_ASSERT(n == 0);

    // Nothing else should have changed
    CPPUNIT_ASSERT(freq.size() == 52);
    CPPUNIT_ASSERT(sqFreq == 58);
  }
}

void DataDefsTest::test_forward_backward_sqfreq() {

  vector<datadefs::num_t> data; 
  map<datadefs::num_t,size_t> freq_left;
  size_t sf_left;
  map<datadefs::num_t,size_t> freq_right;
  size_t sf_right; 
  size_t nRealValues;
  datadefs::num_t x_n;
  size_t n_left;
  size_t n_right;
  for (int i = 0; i < 50; ++i) {
    data.push_back(static_cast<datadefs::num_t>(i));
  }
  data.push_back(0.0);
   
  datadefs::sqfreq(data, freq_left, sf_left, nRealValues);
  datadefs::sqfreq(data, freq_right, sf_right, nRealValues);
  CPPUNIT_ASSERT(freq_left.size() == 50); 
  CPPUNIT_ASSERT(sf_left == 53);
  CPPUNIT_ASSERT(freq_right.size() == 50); 
  CPPUNIT_ASSERT(sf_right == 53);
  
  // Ensure NaN repudiation
  //  !! TODO Eventually, this should be replaced with error code validation. 
  x_n = datadefs::NUM_NAN;
   n_left = 50; n_right = 50;
  datadefs::forward_backward_sqfreq(x_n, n_left, freq_left, sf_left, n_right, freq_right, sf_right);
  
  // Nothing should have changed
  CPPUNIT_ASSERT(n_left == 50);
  CPPUNIT_ASSERT(n_right == 50);
  CPPUNIT_ASSERT(freq_left.size() == 50); 
  CPPUNIT_ASSERT(sf_left == 53);
  CPPUNIT_ASSERT(freq_right.size() == 50); 
  CPPUNIT_ASSERT(sf_right == 53);

  x_n = 10; n_left = 50; n_right = 50;
  datadefs::forward_backward_sqfreq(x_n, n_left, freq_left, sf_left, n_right, freq_right, sf_right);

  CPPUNIT_ASSERT(n_left == 51); // Incremented
  CPPUNIT_ASSERT(freq_left.size() == 50);
  CPPUNIT_ASSERT(sf_left == 56);
  CPPUNIT_ASSERT(n_right == 49); // Decremented
  CPPUNIT_ASSERT(freq_right.size() == 49);
  CPPUNIT_ASSERT(sf_right == 52);
  
  try {
   x_n = 3.1415;
    n_left = static_cast<size_t>(-1); n_right = 0; // Note: n_right would
                                                        // underflow if we
                                                        // didn't bomb out on
                                                        // the overflow first.
    datadefs::forward_backward_sqfreq(x_n, n_left, freq_left, sf_left, n_right, freq_right, sf_right);
    CPPUNIT_FAIL("datadefs::forward_backward_sqfreq didn't throw any exception; expected 'ERRNO_NUMERIC_OVERFLOW'");
  } catch (int e) {  // Ensure standard error is zeroed in the case of an
                     // unusable value for n
    CPPUNIT_ASSERT(e == ERRNO_NUMERIC_OVERFLOW);
    CPPUNIT_ASSERT(n_left == 0); // Overflow

    // Nothing else should have changed
    CPPUNIT_ASSERT(freq_left.size() == 50);
    CPPUNIT_ASSERT(sf_left == 56);
    CPPUNIT_ASSERT(n_right == 0);
    CPPUNIT_ASSERT(freq_right.size() == 49);
    CPPUNIT_ASSERT(sf_right == 52);
  }

  try {
   x_n = 3.1415;
    n_left = 1; n_right = 0; 
    datadefs::forward_backward_sqfreq(x_n, n_left, freq_left, sf_left, n_right, freq_right, sf_right);
    CPPUNIT_FAIL("datadefs::forward_backward_sqfreq didn't throw any exception; expected 'ERRNO_NUMERIC_UNDERFLOW'");
  } catch (int e) {  // Ensure standard error is zeroed in the case of an
                     // unusable value for n
    CPPUNIT_ASSERT(e == ERRNO_NUMERIC_UNDERFLOW);
    CPPUNIT_ASSERT(n_right == static_cast<size_t>(-1)); // Underflow
    CPPUNIT_ASSERT(n_left == 2); // Incremented

    // Nothing else should have changed
    CPPUNIT_ASSERT(freq_left.size() == 50);
    CPPUNIT_ASSERT(sf_left == 56);
    CPPUNIT_ASSERT(freq_right.size() == 49);
    CPPUNIT_ASSERT(sf_right == 52);
  }
}

void DataDefsTest::test_range() {
  vector<size_t> data(50,0);

  datadefs::range(data);
  for (size_t i = 0; i < 50; ++i) {
    CPPUNIT_ASSERT(data[i] == i);
  }

  data.clear();
  datadefs::range(data);
  CPPUNIT_ASSERT(data.size() == 0);

  vector<size_t> atad(50,datadefs::NUM_NAN);
  datadefs::range(atad);
  for (size_t i = 0; i < 50; ++i) {
    CPPUNIT_ASSERT(atad[i] == i);
  }
}

void DataDefsTest::test_sortDataAndMakeRef() {
  vector<datadefs::num_t> data;
  vector<size_t> refIcs;
  data.push_back(0.0);
  for (int i = 0; i < 50; ++i) {
    data.push_back(static_cast<datadefs::num_t>(i));
  }
  
  for (int i = 49; i > -1; --i) { // Deliberately of length data.size() - 1
    refIcs.push_back(static_cast<datadefs::num_t>(i));
  } // This allocation should be irrelevant. We keep it here deliberately to
    //  ensure the results are flattened in a safe manner; we expect a bounds
    //  checker to complain violently should that not be the case. 

  datadefs::sortDataAndMakeRef(true, data, refIcs);
  CPPUNIT_ASSERT(data[0] == 0.0);
  CPPUNIT_ASSERT(refIcs[0] == 0);
  for (int i = 1; i < 51; ++i) {
    CPPUNIT_ASSERT(data[i] == static_cast<datadefs::num_t>(i-1));
    CPPUNIT_ASSERT(refIcs[i] == i);
  }

  datadefs::sortDataAndMakeRef(false, data, refIcs);
  for (int i = 0; i < 50; ++i) {
    CPPUNIT_ASSERT(data[i] == static_cast<datadefs::num_t>(49-i));
    CPPUNIT_ASSERT(refIcs[i] == 50-i);
  }
  CPPUNIT_ASSERT(data[50] == 0.0);
  CPPUNIT_ASSERT(refIcs[50] == 0);

  // Check for correct behavior with an empty data list and arbitrary refIcs
  data.clear();
  datadefs::sortDataAndMakeRef(true, data, refIcs);
  CPPUNIT_ASSERT(data.size() == 0);
  CPPUNIT_ASSERT(refIcs.size() == 0);
  
  datadefs::sortDataAndMakeRef(false, data, refIcs);
  CPPUNIT_ASSERT(data.size() == 0);
  CPPUNIT_ASSERT(refIcs.size() == 0);

  // NaNs are not checked as sorting targets, as their behavior is currently undefined
}

void DataDefsTest::test_utest() {
  vector<datadefs::num_t> x;
  vector<datadefs::num_t> y;
  vector<datadefs::num_t> y2;
  datadefs::num_t pvalue;
  
  for (int i = 0; i < 50; ++i) {
    x.push_back(static_cast<datadefs::num_t>(i));
    y.push_back(static_cast<datadefs::num_t>(i));
    y2.push_back(static_cast<datadefs::num_t>(49-i));
  }

  datadefs::utest(x, y, pvalue);
  CPPUNIT_ASSERT(pvalue == 1);

  datadefs::utest(x, y2, pvalue);
  CPPUNIT_ASSERT(pvalue == 1);

  // Interleave the original input with NaNs; verify we get different results
  for (int i = 0; i < 50; ++i) {
    x.insert(x.begin() + (i*2), datadefs::NUM_NAN);
    y.insert(y.begin() + (i*2), datadefs::NUM_NAN);
    y2.insert(y2.begin() + (i*2), datadefs::NUM_NAN);
  }
  
  datadefs::utest(x, y, pvalue);
  CPPUNIT_ASSERT(fabs(pvalue - 1)
                  < datadefs::EPS);

  datadefs::utest(x, y2, pvalue);
  CPPUNIT_ASSERT(fabs(pvalue - 1)
                  < datadefs::EPS);

  // Verify behavior when x and y are empty
  x.clear();
  y.clear();
  y2.clear();
  
  datadefs::utest(x, y, pvalue);
  CPPUNIT_ASSERT(datadefs::isNAN(pvalue));

  // Verify behavior with all values as NaN
  for (int i = 0; i < 50; ++i) {
    x.push_back(datadefs::NUM_NAN);
    y.push_back(datadefs::NUM_NAN);
  }

  datadefs::utest(x, y, pvalue);
  CPPUNIT_ASSERT(datadefs::isNAN(pvalue));

  // Verify behavior with all values for y as NaN
  x.clear();
  for (int i = 0; i < 50; ++i) {
    x.push_back(static_cast<datadefs::num_t>(i));
  }

  datadefs::utest(x, y, pvalue);
  CPPUNIT_ASSERT(datadefs::isNAN(pvalue));
}

void DataDefsTest::test_erf() {
  datadefs::num_t x;

  x = 0.0;
  CPPUNIT_ASSERT(datadefs::erf(x) == 0.0);

  x = 1.0;
  CPPUNIT_ASSERT(fabs(datadefs::erf(x) - 0.84292558224862812465971728670410811901092529296875)
                  < datadefs::EPS);

  x = datadefs::NUM_NAN;
  CPPUNIT_ASSERT(datadefs::isNAN(datadefs::erf(x)));
}

void DataDefsTest::test_spearman_correlation() {
  CPPUNIT_FAIL("+ This test for datadefs::spearman_correlation(...) is currently unimplemented");
}

void DataDefsTest::test_pearson_correlation() {
  vector<datadefs::num_t> x;
  vector<datadefs::num_t> y;
  vector<datadefs::num_t> y2;
  datadefs::num_t corr;

  for (int i = 0; i < 50; ++i) {
    x.push_back(static_cast<datadefs::num_t>(i));
    y.push_back(static_cast<datadefs::num_t>(i));
    y2.push_back(static_cast<datadefs::num_t>(49-i));
  }
  datadefs::pearson_correlation(x, y, corr);
  CPPUNIT_ASSERT(corr == 1.0);

  datadefs::pearson_correlation(x, y2, corr);
  CPPUNIT_ASSERT(corr == -1.0);
  
  // x and y are defined as const in our signature. Since we're not
  //  testing edge cases of non-trivial memory corruption, we ignore it here.

  x.clear();
  y.clear();
  y2.clear();
  
  datadefs::pearson_correlation(x, y, corr);
  CPPUNIT_ASSERT(datadefs::isNAN(corr));

  for (int i = 0; i < 50; ++i) {
    x.push_back(datadefs::NUM_NAN);
    y.push_back(datadefs::NUM_NAN);
  }
  datadefs::pearson_correlation(x, y, corr);
  CPPUNIT_ASSERT(datadefs::isNAN(corr));
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
                   datadefs::NUM_NAN));
  
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
  justANaN.push_back(datadefs::NUM_NAN);
  CPPUNIT_ASSERT(datadefs::containsNAN(justANaN));

  vector<datadefs::num_t> nanAtBack;
  nanAtBack.push_back(0.0);
  nanAtBack.push_back(1.0);
  nanAtBack.push_back(datadefs::NUM_NAN);
  CPPUNIT_ASSERT(datadefs::containsNAN(nanAtBack));

  vector<datadefs::num_t> nanAtFront;
  nanAtFront.push_back(datadefs::NUM_NAN);
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
  x_n = datadefs::NUM_NAN; n = 0; mu = 0.0; se = 0.0;
  datadefs::forward_sqerr(x_n, n, mu, se);

  // Nothing should have changed
  CPPUNIT_ASSERT(n == 0);
  CPPUNIT_ASSERT(mu == 0.0); 
  CPPUNIT_ASSERT(se == 0.0);

  // Ensure standard error is zeroed in the case of an unusable value for n
  x_n = 3.1415; n = 0; mu = 0.1337; se = 0.7331;
  datadefs::forward_sqerr(x_n, n, mu, se);

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
    CPPUNIT_ASSERT(n == 0); // Overflow

    // Nothing else should have changed

    CPPUNIT_ASSERT(mu == 0.1337);
    CPPUNIT_ASSERT(se == 0.7331);
  }

  x_n = 3.1415; n = 2; mu = 0.1337; se = 0.7331;
  datadefs::forward_sqerr(x_n, n, mu, se);

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
  x_n = datadefs::NUM_NAN;
   n_left = 0; mu_left = 0.0; se_left = 0.0;
   n_right = 0; mu_right = 0.0; se_right = 0.0;
   
  datadefs::forward_backward_sqerr(x_n, n_left, mu_left, se_left, n_right, mu_right, se_right);

  // Nothing should have changed
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
    CPPUNIT_ASSERT(n_left == 0); // Overflow

    // Nothing else should have changed
    CPPUNIT_ASSERT(mu_left == 0.1337);
    CPPUNIT_ASSERT(se_left == 0.7331);
    CPPUNIT_ASSERT(n_right == 0);
    CPPUNIT_ASSERT(mu_right == 0.9773);
    CPPUNIT_ASSERT(se_right == 0.3779);
  }

  try {
   x_n = 3.1415;
    n_left = 1; mu_left = 0.1337; se_left = 0.7331;
    n_right = 0; mu_right = 0.9773; se_right = 0.3779;
    datadefs::forward_backward_sqerr(x_n, n_left, mu_left, se_left, n_right, mu_right, se_right);
    CPPUNIT_FAIL("datadefs::forward_backward_sqerr didn't throw any exception; expected 'ERRNO_NUMERIC_UNDERFLOW'");
  } catch (int e) {  // Ensure standard error is zeroed in the case of an
                     // unusable value for n
    CPPUNIT_ASSERT(e == ERRNO_NUMERIC_UNDERFLOW);
    CPPUNIT_ASSERT(n_right == static_cast<size_t>(-1)); // Underflow
    CPPUNIT_ASSERT(n_left == 2); // Incremented

    // Nothing else should have changed
    CPPUNIT_ASSERT(mu_left == 0.1337);
    CPPUNIT_ASSERT(se_left == 0.7331);
    CPPUNIT_ASSERT(mu_right == 0.9773);
    CPPUNIT_ASSERT(se_right == 0.3779);
  }

  x_n = 3.1415;
   n_left = 2; mu_left = 0.1337; se_left = 0.7331;
   n_right = 4; mu_right = 0.9773; se_right = 0.3779; 

  datadefs::forward_backward_sqerr(x_n, n_left, mu_left, se_left, n_right, mu_right, se_right);

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
                               datadefs::NUM_NAN,0));

  sort(pairedTestVector.begin(), pairedTestVector.end(), incOrder);
  CPPUNIT_ASSERT(pairedTestVector[0].first == -1.0);
  CPPUNIT_ASSERT(pairedTestVector[0].second == 10);
  for (int i = 1; i < 7; ++i) {
    CPPUNIT_ASSERT(pairedTestVector[i].first == static_cast<double>(i-1));
    CPPUNIT_ASSERT(pairedTestVector[i].second == i);
  }
  CPPUNIT_ASSERT(datadefs::isNAN(pairedTestVector[7].first));
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
                               datadefs::NUM_NAN,0));

  sort(pairedTestVector.begin(), pairedTestVector.end(), decOrder);
  for (int i = 0; i > 6; --i) {
    CPPUNIT_ASSERT(pairedTestVector[i].first == static_cast<double>(5-i));
    CPPUNIT_ASSERT(pairedTestVector[i].second == 5-i);
  }
  CPPUNIT_ASSERT(pairedTestVector[6].first == -1.0);
  CPPUNIT_ASSERT(pairedTestVector[6].second == 10);
  CPPUNIT_ASSERT(datadefs::isNAN(pairedTestVector[7].first));
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
  pairedTestVector.push_back(pair<datadefs::num_t, size_t>(datadefs::NUM_NAN,0));

  sort(pairedTestVector.begin(), pairedTestVector.end(), freqIncOrder);
  CPPUNIT_ASSERT(datadefs::isNAN(pairedTestVector[0].first));
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

void DataDefsTest::test_ttest() {
  CPPUNIT_FAIL("+ This test for (deprecated) datadefs::ttest(...) is currently unimplemented");
}

void DataDefsTest::test_regularized_betainc() {
  CPPUNIT_FAIL("+ This test for (deprecated) datadefs::regularized_betainc(...) is currently unimplemented");
}

void DataDefsTest::test_percentile() {
  CPPUNIT_FAIL("+ This test for (deprecated) datadefs::percentile(...) is currently unimplemented");
}

// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( DataDefsTest );

#endif // DATADEFSTEST_HPP
