#ifndef UTILSTEST_HPP
#define UTILSTEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "utils.hpp"

class UtilsTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( UtilsTest );
  CPPUNIT_TEST( test_removeNANs );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();
  
  void test_removeNANs();

};

void UtilsTest::setUp() {}
void UtilsTest::tearDown() {}

void UtilsTest::test_removeNANs() {
  
}

// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( UtilsTest );

#endif // UTILSTEST_HPP
