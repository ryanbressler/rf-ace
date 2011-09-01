#ifndef NODETEST_HPP
#define NODETEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "datadefs.hpp"
#include "node.hpp"
#include "errno.hpp"

class NodeTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE();
  CPPUNIT_TEST( test_setSplitter );
  CPPUNIT_TEST( test_setSplitter );
  CPPUNIT_TEST( test_getSplitter );
  CPPUNIT_TEST( test_percolateData );
  CPPUNIT_TEST( test_getLeafTrainPrediction );
  CPPUNIT_TEST( test_hasChildren );
  CPPUNIT_TEST( test_leafMean );
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
  void setUp();
  void tearDown();

  void test_setSplitter();
  void test_setSplitter();
  void test_getSplitter();
  void test_percolateData();
  void test_getLeafTrainPrediction();
  void test_hasChildren();
  void test_leafMean();
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
  /*
  vector<datadefs::num_t> data;
  datadefs::num_t mu;
  size_t nRealValues;

  for (int i = 0; i < 50; ++i) {
    data.push_back(static_cast<datadefs::num_t>(i));
  }

  // Spuriously assign the values of mu and nRealValues, since they'll be
  //  flattened during function invocation
  mu = -1.0;
  nRealValues = static_cast<size_t>(-1);

  Node::leafMean(data, mu, nRealValues);
  CPPUNIT_ASSERT(mu == 24.5);
  CPPUNIT_ASSERT(nRealValues == 50);

  // Our data vector is defined as const in our signature. Since we're not
  //  testing edge cases of non-trivial memory corruption, we ignore it here.

  // Interleave the original input with NaNs; verify we get the same results
  for (int i = 0; i < 50; ++i) {
    data.insert(data.begin() + (i*2), datadefs::NUM_NAN);
  }

  mu = -1.0;
  nRealValues = static_cast<size_t>(-1);

  Node::leafMean(data, mu, nRealValues);
  CPPUNIT_ASSERT(mu == 24.5);
  CPPUNIT_ASSERT(nRealValues == 50);

  // Ensure a vector containing only NaNs is handled properly
  data.clear();
  for (int i = 0; i < 50; ++i) {
    data.push_back(datadefs::NUM_NAN);
  }

  mu = -1.0;
  nRealValues = static_cast<size_t>(-1);

  Node::leafMean(data, mu, nRealValues);
  CPPUNIT_ASSERT(mu == 0.0);
  CPPUNIT_ASSERT(nRealValues == 0);
  */
}

void NodeTest::test_leafMode() {
  /*
  vector<datadefs::num_t> data;
  datadefs::num_t leafMode = -1.0;
  size_t numClasses = static_cast<size_t>(-1);

  data.push_back(static_cast<datadefs::num_t>(10));
  for (int i = 0; i < 50; ++i) {
    data.push_back(static_cast<datadefs::num_t>(i));
  }

  Node::leafMode(data, leafMode, numClasses);
  CPPUNIT_ASSERT(leafMode == 10.0);

  // numClasses is defined as const in our signature. Since we're not testing
  //  edge cases of non-trivial memory corruption, we ignore it here.

  // Interleave the original input with NaNs; verify we get the same results
  for (int i = 0; i < 50; ++i) {
    data.insert(data.begin() + (i*2), datadefs::NUM_NAN);
  }

  leafMode = -1.0;
  numClasses = static_cast<size_t>(-1);

  Node::leafMode(data, leafMode, numClasses);
  CPPUNIT_ASSERT(leafMode == 10.0);

  // Ensure a vector containing only NaNs is handled as expected
  data.clear();
  for (int i = 0; i < 50; ++i) {
    data.push_back(datadefs::NUM_NAN);
  }

  leafMode = -1.0;
  numClasses = static_cast<size_t>(-1);

  Node::leafMode(data, leafMode, numClasses);
  CPPUNIT_ASSERT(leafMode == 0.0);
  */
}

void NodeTest::test_leafGamma() {
  /*
  vector<datadefs::num_t> data;
  datadefs::num_t leafGamma = -1.0;
  size_t numClasses = 5;

  for (int i = 0; i < 50; ++i) {
    data.push_back(static_cast<datadefs::num_t>(i));
  }

  Node::leafGamma(data, leafGamma, numClasses);
  CPPUNIT_ASSERT(leafGamma == -0.025);

  // numClasses is defined as const in our signature. Since we're not testing
  //  edge cases of non-trivial memory corruption, we ignore it here.

  // Interleave the original input with NaNs; verify we get the same results
  for (int i = 0; i < 50; ++i) {
    data.insert(data.begin() + (i*2), datadefs::NUM_NAN);
  }

  leafGamma = -1.0;
  numClasses = 5;

  Node::leafGamma(data, leafGamma, numClasses);
  CPPUNIT_ASSERT(leafGamma == -0.025);

  // Ensure a vector containing only NaNs is handled as expected
  data.clear();
  for (int i = 0; i < 50; ++i) {
    data.push_back(datadefs::NUM_NAN);
  }

  leafGamma = -1.0;
  numClasses = 5;

  Node::leafGamma(data, leafGamma, numClasses);
  CPPUNIT_ASSERT(leafGamma == 0.0);
  */
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
