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
  //CPPUNIT_TEST( test_strv2catv );
  //CPPUNIT_TEST( test_strv2numv );
  CPPUNIT_TEST( test_cardinality );
  CPPUNIT_TEST( test_countRealValues );
  CPPUNIT_TEST( test_map_data );
  //CPPUNIT_TEST( test_max_idx );
  //CPPUNIT_TEST( test_sortDataAndMakeRef );
  CPPUNIT_TEST( test_isNAN );
  CPPUNIT_TEST( test_containsNAN );
  CPPUNIT_TEST( test_increasingOrderOperator );
  CPPUNIT_TEST( test_decreasingOrderOperator );
  CPPUNIT_TEST( test_freqIncreasingOrderOperator );
  CPPUNIT_TEST( test_make_pairedv );
  CPPUNIT_TEST( test_separate_pairedv );
  //CPPUNIT_TEST( test_sortFromRef );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();
  
  void test_size_tIsSigned();
  //void test_strv2catv();
  //void test_strv2numv();
  void test_cardinality();
  void test_countRealValues();
  void test_map_data();
  //void test_max_idx();
  //void test_sortDataAndMakeRef();
  void test_percentile();
  void test_isNAN();
  void test_containsNAN();
  void test_increasingOrderOperator();
  void test_decreasingOrderOperator();
  void test_freqIncreasingOrderOperator();
  void test_make_pairedv();
  void test_separate_pairedv();
  //void test_sortFromRef();
};

void DataDefsTest::setUp() {}
void DataDefsTest::tearDown() {}

void DataDefsTest::test_size_tIsSigned() {
  
  // Ensure size_t is conformant with the 1999 ISO C standard (C99)
  CPPUNIT_ASSERT(!numeric_limits<size_t>::is_signed);
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

void DataDefsTest::test_isNAN() {

  // Test for isNAN(&string)
  CPPUNIT_ASSERT(datadefs::isNAN_STR("NA"));
  CPPUNIT_ASSERT(datadefs::isNAN_STR("NAN"));
  CPPUNIT_ASSERT(datadefs::isNAN_STR("?"));

  CPPUNIT_ASSERT(!datadefs::isNAN_STR("2"));
  CPPUNIT_ASSERT(!datadefs::isNAN_STR("@data"));
  CPPUNIT_ASSERT(!datadefs::isNAN_STR("NAte"));

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

// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( DataDefsTest );

#endif // DATADEFSTEST_HPP
