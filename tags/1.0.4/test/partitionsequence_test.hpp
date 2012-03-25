#ifndef PARTITIONSEQUENCETEST_HPP
#define PARTITIONSEQUENCETEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "partitionsequence.hpp"
#include "errno.hpp"

class PartitionSequenceTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( PartitionSequenceTest );
  CPPUNIT_TEST( test_partitionSequenceOfLength4 );
  CPPUNIT_TEST( test_partitionSequenceOfLength16 );
  //CPPUNIT_TEST( test_isNextBitAdded );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();
  void test_partitionSequenceOfLength4();
  void test_partitionSequenceOfLength16();
  //void test_isNextBitAdded();

};

void PartitionSequenceTest::setUp() {}
void PartitionSequenceTest::tearDown() {}

void PartitionSequenceTest::test_partitionSequenceOfLength4() {
  
  size_t sequenceLength = 4;
  
  PartitionSequence::PartitionSequence PS;
  
  CPPUNIT_ASSERT( PS.at(0) == 0 && PS.isAdded(0)  == true);
  CPPUNIT_ASSERT( PS.at(1) == 1 && PS.isAdded(1)  == true);
  CPPUNIT_ASSERT( PS.at(2) == 0 && PS.isAdded(2)  == false);
  CPPUNIT_ASSERT( PS.at(3) == 2 && PS.isAdded(3)  == true);
  CPPUNIT_ASSERT( PS.at(4) == 0 && PS.isAdded(4)  == true);
  CPPUNIT_ASSERT( PS.at(5) == 1 && PS.isAdded(5)  == false);
  CPPUNIT_ASSERT( PS.at(6) == 0 && PS.isAdded(6)  == false);
  CPPUNIT_ASSERT( PS.at(7) == 3 && PS.isAdded(7)  == true);
  CPPUNIT_ASSERT( PS.at(8) == 0 && PS.isAdded(8)  == true);
  CPPUNIT_ASSERT( PS.at(9) == 1 && PS.isAdded(9)  == true);
  CPPUNIT_ASSERT( PS.at(10) == 0 && PS.isAdded(10) == false);
  CPPUNIT_ASSERT( PS.at(11) == 2 && PS.isAdded(11) == false);
  CPPUNIT_ASSERT( PS.at(12) == 0 && PS.isAdded(12) == true);
  CPPUNIT_ASSERT( PS.at(13) == 1 && PS.isAdded(13) == false);
  CPPUNIT_ASSERT( PS.at(14) == 0 && PS.isAdded(14) == false);
  
}

void PartitionSequenceTest::test_partitionSequenceOfLength16() {

  size_t sequenceLength = 40;

  PartitionSequence::PartitionSequence PS;

  // 2^( 16 - 2) - 1
  size_t psMax = ( 1 << ( 16 - 2 ) ) - 1;

  //psMax = 2^16 - 1;

  for ( size_t psIdx = 0; psIdx <= psMax; ++psIdx ) {
    CPPUNIT_ASSERT( PS.at(psIdx) < 16 );
  }

}

// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( PartitionSequenceTest ); 

#endif // PARTITIONSEQUENCETEST_HPP
