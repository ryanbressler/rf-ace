#ifndef ROOTNODETEST_HPP
#define ROOTNODETEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "datadefs.hpp"
#include "rootnode.hpp"
#include "errno.hpp"

class RootNodeTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( RootNodeTest );
  CPPUNIT_TEST( test_getTreeSizeEstimate );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();
  void test_getTreeSizeEstimate();

};

void RootNodeTest::setUp() {}
void RootNodeTest::tearDown() {}

void RootNodeTest::test_getTreeSizeEstimate() {

  RootNode rootNode;

  size_t nSamples;
  size_t nMaxLeaves = 1;
  size_t nodeSize = 1;

  for ( size_t nSamples = 1; nSamples < 10; ++nSamples ) {
    CPPUNIT_ASSERT( rootNode.getTreeSizeEstimate(nSamples,nMaxLeaves,nodeSize) == 1 );
  }

}

// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( RootNodeTest ); 

#endif // ROOTNODETEST_HPP
