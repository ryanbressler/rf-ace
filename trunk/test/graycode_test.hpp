#ifndef GRAYCODETEST_HPP
#define GRAYCODETEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "graycode.hpp"
#include "errno.hpp"

class GrayCodeTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( GrayCodeTest );
  //CPPUNIT_TEST( test_setSplitter );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();

};

void GrayCodeTest::setUp() {}
void GrayCodeTest::tearDown() {}

// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( GrayCodeTest ); 

#endif // GRAYCODETEST_HPP
