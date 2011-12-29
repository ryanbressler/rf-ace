#ifndef SPLITTERTEST_HPP
#define SPLITTERTEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "splitter.hpp"
#include "errno.hpp"

class splitterTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( splitterTest );
  CPPUNIT_TEST( test_NumericalSplitterSplitsLeftAndRight );
  CPPUNIT_TEST( test_CategoricalSplitterSplitsLeftAndRight );
  //CPPUNIT_TEST( test_isNextBitAdded );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();
  void test_NumericalSplitterSplitsLeftAndRight();
  void test_CategoricalSplitterSplitsLeftAndRight();
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

  
  set<num_t> splitValuesLeft;
  splitValuesLeft.insert(0);
  splitValuesLeft.insert(1);
  splitValuesLeft.insert(3);
  
  // NOTE: value 4 is missing from the union of left and right split sets
  
  set<num_t> splitValuesRight;
  splitValuesRight.insert(2);
  splitValuesRight.insert(5);
  
  Splitter::Splitter splitter("foo",splitValuesLeft,splitValuesRight);
  
  CPPUNIT_ASSERT( splitter.splitterName_ == "foo" );
  CPPUNIT_ASSERT( splitter.name() == "foo" );
    
  CPPUNIT_ASSERT( splitter.splitsLeft(0) );
  CPPUNIT_ASSERT( splitter.splitsLeft(1) );
  CPPUNIT_ASSERT( splitter.splitsLeft(3) );
  
  CPPUNIT_ASSERT( !splitter.splitsLeft(2) );
  CPPUNIT_ASSERT( !splitter.splitsLeft(5) );
  
  CPPUNIT_ASSERT( !splitter.splitsRight(0) );
  CPPUNIT_ASSERT( !splitter.splitsRight(1) );
  CPPUNIT_ASSERT( !splitter.splitsRight(3) );
  
  CPPUNIT_ASSERT( splitter.splitsRight(2) );
  CPPUNIT_ASSERT( splitter.splitsRight(5) );
    
  CPPUNIT_ASSERT( !splitter.splitsLeft(4) );
  CPPUNIT_ASSERT( !splitter.splitsRight(4) );
  

  
  
}


// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( splitterTest ); 

#endif // SPLITTERTEST_HPP
