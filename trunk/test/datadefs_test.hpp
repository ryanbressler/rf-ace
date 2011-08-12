#ifndef DATADEFSTEST_HPP
#define DATADEFSTEST_HPP

#include "argparse.hpp"
#include "errno.hpp"
#include <cppunit/extensions/HelperMacros.h>

class DataDefsTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( DataDefsTest );
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

void DataDefsTest::test_strv2catv() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_strv2numv() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_str2num() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_meanVals() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_mean() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_mode() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_gamma() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_cardinality() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_sqerr() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_countRealValues() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_count_freq() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_map_data() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_gini() {
  // !! Note: two signatures!
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_sqfreq() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_forward_sqfreq() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_forward_backward_sqfreq() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_range() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_sortDataAndMakeRef() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_ttest() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_utest() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_erf() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_regularized_betainc() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_spearman_correlation() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_pearson_correlation() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_percentile() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_isNAN() {
  // !! Note: two (or more) signatures!
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_containsNAN() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_forward_sqerr() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_forward_backward_sqerr() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_increasingOrderOperator() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_decreasingOrderOperator() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_freqIncreasingOrderOperator() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_make_pairedv() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_separate_pairedv() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

void DataDefsTest::test_sortFromRef() {
  CPPUNIT_FAIL("+ This test is currently unimplemented");
}

// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( DataDefsTest );

#endif // DATADEFSTEST_HPP
