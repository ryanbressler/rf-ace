#ifndef SPLITTERTEST_HPP
#define SPLITTERTEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "splitter.hpp"
#include "treedata.hpp"
#include "errno.hpp"

class splitterTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( splitterTest );
  CPPUNIT_TEST( test_NumericalSplitterSplitsLeftAndRight );
  CPPUNIT_TEST( test_CategoricalSplitterSplitsLeftAndRight );
  CPPUNIT_TEST( test_NumericalSplittingWithNovelData );
  CPPUNIT_TEST( test_CategoricalSplittingWithNovelData );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();
  void test_NumericalSplitterSplitsLeftAndRight();
  void test_CategoricalSplitterSplitsLeftAndRight();
  void test_NumericalSplittingWithNovelData();
  void test_CategoricalSplittingWithNovelData();
};

void splitterTest::setUp() {}
void splitterTest::tearDown() {}
 
void splitterTest::test_NumericalSplitterSplitsLeftAndRight() {

  num_t splitLeftLeqValue = 0.5;

  Splitter::Splitter splitter("foo",splitLeftLeqValue);

  CPPUNIT_ASSERT( splitter.splitterName_ == "foo" );
  CPPUNIT_ASSERT( splitter.name() == "foo" );

  CPPUNIT_ASSERT( splitter.splitsLeft(0.3) );
  CPPUNIT_ASSERT( !splitter.splitsLeft(0.66) );

  CPPUNIT_ASSERT( splitter.splitsRight(0.99) );
  CPPUNIT_ASSERT( !splitter.splitsRight(-0.1) );
  
}

void splitterTest::test_CategoricalSplitterSplitsLeftAndRight() {

  
  set<string> splitValuesLeft;
  splitValuesLeft.insert("0");
  splitValuesLeft.insert("1");
  splitValuesLeft.insert("3");
  
  // NOTE: value 4 is missing from the union of left and right split sets
  
  set<string> splitValuesRight;
  splitValuesRight.insert("2");
  splitValuesRight.insert("5");
  
  Splitter::Splitter splitter("foo",splitValuesLeft,splitValuesRight);
  
  CPPUNIT_ASSERT( splitter.splitterName_ == "foo" );
  CPPUNIT_ASSERT( splitter.name() == "foo" );
    
  CPPUNIT_ASSERT( splitter.splitsLeft("0") );
  CPPUNIT_ASSERT( splitter.splitsLeft("1") );
  CPPUNIT_ASSERT( splitter.splitsLeft("3") );
  
  CPPUNIT_ASSERT( !splitter.splitsLeft("2") );
  CPPUNIT_ASSERT( !splitter.splitsLeft("5") );
  
  CPPUNIT_ASSERT( !splitter.splitsRight("0") );
  CPPUNIT_ASSERT( !splitter.splitsRight("1") );
  CPPUNIT_ASSERT( !splitter.splitsRight("3") );
  
  CPPUNIT_ASSERT( splitter.splitsRight("2") );
  CPPUNIT_ASSERT( splitter.splitsRight("5") );
    
  CPPUNIT_ASSERT( !splitter.splitsLeft("4") );
  CPPUNIT_ASSERT( !splitter.splitsRight("4") );
  
}

void splitterTest::test_NumericalSplittingWithNovelData() {
  
  
  Treedata treeData("test_6by10_mixed_matrix.tsv",'\t',':');

  /*
    foo     S1      S2      S3      S4      S5      S6      S7      S8      S9      S10
    N:F1    nA      8.5     3.4     7.2     5       6       7       11      9       NA
    N:F2    2       3       4       5       6       NA      NA      9       nan     10
    C:F3    NA      nA      naN     NaN     1       1       1       2       2       2
    N:F4    10      9.9     8       7       6       5       4       3       2.4     1
    C:F5    3       3       3       4       4       5       3       2       2       2
    N:F6    9       8       7       9       8       7       3       2       1.0     99.23
  */    

  Splitter splitter1("N:F1",5.1);

  // samples 0 and 9 (NA) Should not split at all
  CPPUNIT_ASSERT( !splitter1.splitsLeft(&treeData,0) );
  CPPUNIT_ASSERT( !splitter1.splitsRight(&treeData,0) );

  CPPUNIT_ASSERT( !splitter1.splitsLeft(&treeData,9) );
  CPPUNIT_ASSERT( !splitter1.splitsRight(&treeData,9) );

  // splits right: samples 1 3 5 6 7 8 
  CPPUNIT_ASSERT( splitter1.splitsRight(&treeData,1) );
  CPPUNIT_ASSERT( splitter1.splitsRight(&treeData,3) );
  CPPUNIT_ASSERT( splitter1.splitsRight(&treeData,5) );
  CPPUNIT_ASSERT( splitter1.splitsRight(&treeData,6) );
  CPPUNIT_ASSERT( splitter1.splitsRight(&treeData,7) );
  CPPUNIT_ASSERT( splitter1.splitsRight(&treeData,8) );
  CPPUNIT_ASSERT( !splitter1.splitsLeft(&treeData,1) );
  CPPUNIT_ASSERT( !splitter1.splitsLeft(&treeData,3) );
  CPPUNIT_ASSERT( !splitter1.splitsLeft(&treeData,5) );
  CPPUNIT_ASSERT( !splitter1.splitsLeft(&treeData,6) );
  CPPUNIT_ASSERT( !splitter1.splitsLeft(&treeData,7) );
  CPPUNIT_ASSERT( !splitter1.splitsLeft(&treeData,8) );

  // splits left: samples 2 4
  CPPUNIT_ASSERT( splitter1.splitsLeft(&treeData,2) );
  CPPUNIT_ASSERT( splitter1.splitsLeft(&treeData,4) );
  CPPUNIT_ASSERT( !splitter1.splitsRight(&treeData,2) );
  CPPUNIT_ASSERT( !splitter1.splitsRight(&treeData,4) );

}

void splitterTest::test_CategoricalSplittingWithNovelData() {

  Treedata treeData("test_6by10_mixed_matrix.tsv",'\t',':');

  /*
    foo     S1      S2      S3      S4      S5      S6      S7      S8      S9      S10
    N:F1    nA      8.5     3.4     7.2     5       6       7       11      9       NA
    N:F2    2       3       4       5       6       NA      NA      9       nan     10
    C:F3    NA      nA      naN     NaN     1       1       1       2       2       2
    N:F4    10      9.9     8       7       6       5       4       3       2.4     1
    C:F5    3       3       3       4       4       5       3       2       2       2
    N:F6    9       8       7       9       8       7       3       2       1.0     99.23
  */

  // Splitter for C:F3

  set<string> splitValues_left,splitValues_right;

  splitValues_left.insert("1");
  splitValues_right.insert("2");

  Splitter splitter1("C:F3",splitValues_left,splitValues_right);

  // splits neither left nor right (NA): samples 0 1 2 3
  CPPUNIT_ASSERT( !splitter1.splitsLeft(&treeData,0) );
  CPPUNIT_ASSERT( !splitter1.splitsRight(&treeData,0) );
  CPPUNIT_ASSERT( !splitter1.splitsLeft(&treeData,1) );
  CPPUNIT_ASSERT( !splitter1.splitsRight(&treeData,1) );
  CPPUNIT_ASSERT( !splitter1.splitsLeft(&treeData,2) );
  CPPUNIT_ASSERT( !splitter1.splitsRight(&treeData,2) );
  CPPUNIT_ASSERT( !splitter1.splitsLeft(&treeData,3) );
  CPPUNIT_ASSERT( !splitter1.splitsRight(&treeData,3) );

  // splits left: 4 5 6
  CPPUNIT_ASSERT( splitter1.splitsLeft(&treeData,4) );
  CPPUNIT_ASSERT( splitter1.splitsLeft(&treeData,5) );
  CPPUNIT_ASSERT( splitter1.splitsLeft(&treeData,6) );
  CPPUNIT_ASSERT( !splitter1.splitsRight(&treeData,4) );
  CPPUNIT_ASSERT( !splitter1.splitsRight(&treeData,5) );
  CPPUNIT_ASSERT( !splitter1.splitsRight(&treeData,6) );

  // splits right: 7 8 9
  CPPUNIT_ASSERT( !splitter1.splitsLeft(&treeData,7) );
  CPPUNIT_ASSERT( !splitter1.splitsLeft(&treeData,8) );
  CPPUNIT_ASSERT( !splitter1.splitsLeft(&treeData,9) );
  CPPUNIT_ASSERT( splitter1.splitsRight(&treeData,7) );
  CPPUNIT_ASSERT( splitter1.splitsRight(&treeData,8) );
  CPPUNIT_ASSERT( splitter1.splitsRight(&treeData,9) );

  // C:F5    3       3       3       4       4       5       3       2       2       2

  splitValues_left.clear();
  splitValues_left.insert("3");

  splitValues_right.clear();
  splitValues_right.insert("4");

  Splitter splitter2("C:F5",splitValues_left,splitValues_right);

  // splits neither left nor right: samples 5 7 8 9
  CPPUNIT_ASSERT( !splitter2.splitsLeft(&treeData,5) );
  CPPUNIT_ASSERT( !splitter2.splitsRight(&treeData,5) );
  CPPUNIT_ASSERT( !splitter2.splitsLeft(&treeData,7) );
  CPPUNIT_ASSERT( !splitter2.splitsRight(&treeData,7) );
  CPPUNIT_ASSERT( !splitter2.splitsLeft(&treeData,8) );
  CPPUNIT_ASSERT( !splitter2.splitsRight(&treeData,8) );
  CPPUNIT_ASSERT( !splitter2.splitsLeft(&treeData,9) );
  CPPUNIT_ASSERT( !splitter2.splitsRight(&treeData,9) );

}

// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( splitterTest ); 

#endif // SPLITTERTEST_HPP
