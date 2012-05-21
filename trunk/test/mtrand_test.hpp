#ifndef MTRANDTEST_HPP
#define MTRANDTEST_HPP

#include <cppunit/extensions/HelperMacros.h>

class MTRandTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( MTRandTest );
  //CPPUNIT_TEST( test_offByOneSeeds );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();
  //void test_offByOneSeeds();

};

void MTRandTest::setUp() {}
void MTRandTest::tearDown() {}

// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( MTRandTest ); 

#endif // MTRANDTEST_HPP
