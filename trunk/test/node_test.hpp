#ifndef NODETEST_HPP
#define NODETEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "datadefs.hpp"
#include "treedata.hpp"
#include "node.hpp"
#include "errno.hpp"

class NodeTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( NodeTest );
  CPPUNIT_TEST( test_setSplitter );
  CPPUNIT_TEST( test_percolateData );
  CPPUNIT_TEST( test_getLeafTrainPrediction );
  CPPUNIT_TEST( test_hasChildren );
  CPPUNIT_TEST( test_recursiveNodeSplit );
  CPPUNIT_TEST( test_cleanPairVectorFromNANs );
  CPPUNIT_TEST( test_recursiveNDescendantNodes );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();

  void test_setSplitter();
  void test_percolateData();
  void test_getLeafTrainPrediction();
  void test_hasChildren();
  void test_recursiveNodeSplit();
  void test_cleanPairVectorFromNANs();
  void test_recursiveNDescendantNodes();

};

void NodeTest::setUp() {}
void NodeTest::tearDown() {}


void NodeTest::test_setSplitter() {
  
  size_t splitterIdx = 3;
  datadefs::num_t splitLeftLeqValue = 0.5; 
  datadefs::num_t leftFraction = 0.5;
  
  //Splitter::Splitter splitter(0.5);
  
  Node::Node node;
    
  node.setSplitter(splitterIdx,"foo",leftFraction,splitLeftLeqValue);
  
  CPPUNIT_ASSERT( node.splitterIdx() == splitterIdx );
  CPPUNIT_ASSERT( node.splitter_.isNumerical );
  CPPUNIT_ASSERT( fabs(node.splitter_.leftLeqValue - splitLeftLeqValue) < datadefs::EPS );
  CPPUNIT_ASSERT( fabs(node.splitter_.leftFraction - leftFraction) < datadefs::EPS );
  
}

void NodeTest::test_percolateData() {
  
  Node::Node node0;
  //Splitter splitter("foo",0.1);
  node0.setSplitter(1,"foo",0.5,0.1);
  CPPUNIT_ASSERT( node0.leftChild() == node0.percolateData(0.09) );
  CPPUNIT_ASSERT( node0.rightChild() == node0.percolateData(0.11) );
  CPPUNIT_ASSERT( NULL == node0.percolateData(datadefs::NUM_NAN));

  Node::Node node1;
  
  set<string> leftValues;
  set<string> rightValues;

  leftValues.insert("a");
  leftValues.insert("b");

  rightValues.insert("c");
  rightValues.insert("d");

  node1.setSplitter(1,"foo",0.5,leftValues,rightValues);

  CPPUNIT_ASSERT( node1.percolateData("a") == node1.leftChild() );
  CPPUNIT_ASSERT( node1.percolateData("b") == node1.leftChild() );
  CPPUNIT_ASSERT( node1.percolateData("c") == node1.rightChild() );
  CPPUNIT_ASSERT( node1.percolateData("d") == node1.rightChild() );

  CPPUNIT_ASSERT( node1.percolateData("foo") == NULL );

}

void NodeTest::test_getLeafTrainPrediction() {
}

void NodeTest::test_hasChildren() {
}

void NodeTest::test_recursiveNodeSplit() { 
}

void NodeTest::test_cleanPairVectorFromNANs() { 

}

void NodeTest::test_recursiveNDescendantNodes() {
  
}
// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( NodeTest ); 

#endif // NODETEST_HPP
