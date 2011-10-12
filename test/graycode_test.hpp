#ifndef GRAYCODETEST_HPP
#define GRAYCODETEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "graycode.hpp"
#include "errno.hpp"

class GrayCodeTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( GrayCodeTest );
  CPPUNIT_TEST( test_getCode );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();
  void test_getCode();

};

void GrayCodeTest::setUp() {}
void GrayCodeTest::tearDown() {}

void GrayCodeTest::test_getCode() {

  size_t nMaxBits = 4;
  
  GrayCode::GrayCode GC(nMaxBits);
  
  CPPUNIT_ASSERT( GC.getCode(0)  ==  0 );
  CPPUNIT_ASSERT( GC.getCode(1)  ==  1 );
  CPPUNIT_ASSERT( GC.getCode(2)  ==  3 );
  CPPUNIT_ASSERT( GC.getCode(3)  ==  2 );
  CPPUNIT_ASSERT( GC.getCode(4)  ==  6 );
  CPPUNIT_ASSERT( GC.getCode(5)  ==  7 );
  CPPUNIT_ASSERT( GC.getCode(6)  ==  5 );
  CPPUNIT_ASSERT( GC.getCode(7)  ==  4 );
  CPPUNIT_ASSERT( GC.getCode(8)  == 12 );
  CPPUNIT_ASSERT( GC.getCode(9)  == 13 );
  CPPUNIT_ASSERT( GC.getCode(10) == 15 );
  CPPUNIT_ASSERT( GC.getCode(11) == 14 );
  CPPUNIT_ASSERT( GC.getCode(12) == 10 );
  CPPUNIT_ASSERT( GC.getCode(13) == 11 );
  CPPUNIT_ASSERT( GC.getCode(14) ==  9 );
  CPPUNIT_ASSERT( GC.getCode(15) ==  8 );

}

// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( GrayCodeTest ); 

#endif // GRAYCODETEST_HPP
