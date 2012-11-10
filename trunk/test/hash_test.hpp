#ifndef HASH_HPP
#define HASH_HPP

#include <cppunit/extensions/HelperMacros.h>
//#include "datadefs.hpp"
//#include "hash.hpp"
//#include "errno.hpp"

class HashTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( HashTest );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();

};

void HashTest::setUp() {}
void HashTest::tearDown() {}


// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( HashTest ); 

#endif // HASH_HPP
