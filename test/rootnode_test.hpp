#ifndef ROOTNODETEST_HPP
#define ROOTNODETEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "datadefs.hpp"
#include "rootnode.hpp"
#include "errno.hpp"

class RootNodeTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( RootNodeTest );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();

};

void RootNodeTest::setUp() {}
void RootNodeTest::tearDown() {}

// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( RootNodeTest ); 

#endif // ROOTNODETEST_HPP
