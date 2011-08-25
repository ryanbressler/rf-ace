#ifndef NODETEST_HPP
#define NODETEST_HPP

#ifndef TEST__
#define TEST__
#endif

#include "node.hpp"
#include "errno.hpp"
#include <cppunit/extensions/HelperMacros.h>

class NodeTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE();
  CPPUNIT_TEST( test_setSplitter );
  CPPUNIT_TEST( test_setSplitter );
  CPPUNIT_TEST( test_getSplitter );
  CPPUNIT_TEST( test_percolateData );
  CPPUNIT_TEST( test_getLeafTrainPrediction );
  CPPUNIT_TEST( test_hasChildren );  CPPUNIT_TEST( test_leafMean );
  CPPUNIT_TEST( test_leafMode );
  CPPUNIT_TEST( test_leafGamma );
  CPPUNIT_TEST( test_recursiveNodeSplit );
  CPPUNIT_TEST( test_cleanPairVectorFromNANs );
  CPPUNIT_TEST( test_numericalFeatureSplit );
  CPPUNIT_TEST( test_categoricalFeatureSplit );
  CPPUNIT_TEST( test_splitFitness );
  CPPUNIT_TEST( test_recursiveNDescendantNodes );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void test_setUp();
  void test_tearDown();

  void test_setSplitter();
  void test_setSplitter();
  void test_getSplitter();
  void test_percolateData();
  void test_getLeafTrainPrediction();
  void test_hasChildren();  void test_leafMean();
  void test_leafMode();
  void test_leafGamma();
  void test_recursiveNodeSplit();
  void test_cleanPairVectorFromNANs();
  void test_numericalFeatureSplit();
  void test_categoricalFeatureSplit();
  void test_splitFitness();
  void test_recursiveNDescendantNodes();

};

void NodeTest::setUp() {}
void NodeTest::tearDown() {}

void NodeTest::test_setSplitter() {
  // Two signatures!
  //void Node::setSplitter(size_t splitter, num_t threshold);
  //void Node::setSplitter(size_t splitter, set<num_t> classSet);
}

void NodeTest::test_getSplitter() {
  //Node::getSplitter(...);
}

void NodeTest::test_percolateData() {
  //Node::percolateData(...);
}

void NodeTest::test_getLeafTrainPrediction() {
}

void NodeTest::test_hasChildren() {
}

void NodeTest::test_leafMean() {
}

void NodeTest::test_leafMode() { 
}

void NodeTest::test_leafGamma() { 
}

void NodeTest::test_recursiveNodeSplit() { 
}

void NodeTest::test_cleanPairVectorFromNANs() { 
}

void NodeTest::test_numericalFeatureSplit() { 
}

void NodeTest::test_categoricalFeatureSplit() { 
}

void NodeTest::test_splitFitness() { 
}

void NodeTest::test_recursiveNDescendantNodes() {
  
}
// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION() { 

#endif // NODETEST_HPP
