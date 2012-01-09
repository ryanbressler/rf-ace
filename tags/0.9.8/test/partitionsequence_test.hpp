#ifndef PARTITIONSEQUENCETEST_HPP
#define PARTITIONSEQUENCETEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "partitionsequence.hpp"
#include "errno.hpp"

class PartitionSequenceTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( PartitionSequenceTest );
  CPPUNIT_TEST( test_at );
  CPPUNIT_TEST( test_isAdded );
  //CPPUNIT_TEST( test_isNextBitAdded );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();
  void test_at();
  void test_isAdded();
  //void test_isNextBitAdded();

};

void PartitionSequenceTest::setUp() {}
void PartitionSequenceTest::tearDown() {}

void PartitionSequenceTest::test_at() {
  
  size_t sequenceLength = 4;
  
  PartitionSequence::PartitionSequence PS(sequenceLength);
  
  //cout << GC.nextBitDiff(0) << endl;
  
  CPPUNIT_ASSERT( PS.at(0) == 0 );
  CPPUNIT_ASSERT( PS.at(1) == 1 );
  CPPUNIT_ASSERT( PS.at(2) == 0 );
  CPPUNIT_ASSERT( PS.at(3) == 2 );
  CPPUNIT_ASSERT( PS.at(4) == 0 );
  CPPUNIT_ASSERT( PS.at(5) == 1 );
  CPPUNIT_ASSERT( PS.at(6) == 0 );
  CPPUNIT_ASSERT( PS.at(7) == 3 );
  CPPUNIT_ASSERT( PS.at(8) == 0 );
  CPPUNIT_ASSERT( PS.at(9) == 1 );
  CPPUNIT_ASSERT( PS.at(10) == 0 );
  CPPUNIT_ASSERT( PS.at(11) == 2 );
  CPPUNIT_ASSERT( PS.at(12) == 0 );
  CPPUNIT_ASSERT( PS.at(13) == 1 );
  CPPUNIT_ASSERT( PS.at(14) == 0 );
  
}
 
void PartitionSequenceTest::test_isAdded() {
  
  size_t sequenceLength = 4;
  
  PartitionSequence::PartitionSequence PS(sequenceLength);
  
  CPPUNIT_ASSERT( PS.isAdded(0)  == true );
  CPPUNIT_ASSERT( PS.isAdded(1)  == true );
  CPPUNIT_ASSERT( PS.isAdded(2)  == false );
  CPPUNIT_ASSERT( PS.isAdded(3)  == true );
  CPPUNIT_ASSERT( PS.isAdded(4)  == true );
  CPPUNIT_ASSERT( PS.isAdded(5)  == false );
  CPPUNIT_ASSERT( PS.isAdded(6)  == false );
  CPPUNIT_ASSERT( PS.isAdded(7)  == true );
  CPPUNIT_ASSERT( PS.isAdded(8)  == true );
  CPPUNIT_ASSERT( PS.isAdded(9)  == true );
  CPPUNIT_ASSERT( PS.isAdded(10) == false );
  CPPUNIT_ASSERT( PS.isAdded(11) == false );
  CPPUNIT_ASSERT( PS.isAdded(12) == true );
  CPPUNIT_ASSERT( PS.isAdded(13) == false );
  CPPUNIT_ASSERT( PS.isAdded(14) == false );
}

// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( PartitionSequenceTest ); 

#endif // PARTITIONSEQUENCETEST_HPP
