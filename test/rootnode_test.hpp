#ifndef ROOTNODETEST_HPP
#define ROOTNODETEST_HPP

#include <cppunit/extensions/HelperMacros.h>
#include "datadefs.hpp"
#include "partitionsequence.hpp"
#include "rootnode.hpp"
#include "errno.hpp"

class RootNodeTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE( RootNodeTest );
  CPPUNIT_TEST( test_regressionTreeSplitting ); 
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();

  void test_regressionTreeSplitting();

};

void RootNodeTest::setUp() {}
void RootNodeTest::tearDown() {}

void RootNodeTest::test_regressionTreeSplitting() {
    
  bool sampleWithReplacement = true;
  datadefs::num_t sampleSizeFraction = 1.0;
  size_t maxNodesToStop = 100;
  size_t minNodeSizeToStop = 3;
  bool isRandomSplit = true;
  size_t nFeaturesForSplit = 1;
  bool useContrasts = true;
  size_t numClasses = 0;
  PartitionSequence::PartitionSequence* partitionSequence = new PartitionSequence(5);

  RootNode::RootNode rootNode(sampleWithReplacement,
			      sampleSizeFraction,
			      maxNodesToStop,
			      minNodeSizeToStop,
			      isRandomSplit,
			      nFeaturesForSplit,
			      useContrasts,
			      numClasses,
			      partitionSequence);

  // Split the rootnode
  rootNode.setSplitter(0,"foo",1.1);
  rootNode.trainPrediction_ = 5.0;
  Node::Node* node0 = rootNode.leftChild_;
  Node::Node* node1 = rootNode.rightChild_;

  // Set node0 (left child of rootnode)
  set<num_t> splitValuesLeft;
  splitValuesLeft.insert(0);
  splitValuesLeft.insert(1);
  set<num_t> splitValuesRight;
  splitValuesRight.insert(2);
  node0->setSplitter(1,"foo",splitValuesLeft,splitValuesRight);
  node0->trainPrediction_ = 4.0;

  // Set node1 (right child of rootnode)
  node1->setSplitter(2,"foo",3.0);
  node1->trainPrediction_ = 6.0;

  // Children of node0
  Node::Node* node00 = node0->leftChild_;
  Node::Node* node01 = node0->rightChild_; 

  // Set nodes 00 and 01
  node00->trainPrediction_ = 3.9;
  node01->setSplitter(3,"foo",-1.5);
  node01->trainPrediction_ = 4.2;

  // Children of node1
  Node::Node* node10 = node1->leftChild_;
  Node::Node* node11 = node1->rightChild_;

  // Set nodes 10 and 11
  node10->trainPrediction_ = 5.1;
  splitValuesLeft.clear();
  splitValuesLeft.insert(2);
  splitValuesLeft.insert(3);
  splitValuesRight.clear();
  splitValuesRight.insert(0);
  splitValuesRight.insert(1);
  node11->setSplitter(4,"foo",splitValuesLeft,splitValuesRight);
  node11->trainPrediction_ = 6.6;

  // Children of node01
  Node::Node* node010 = node01->leftChild_;
  Node::Node* node011 = node01->rightChild_;

  // Set nodes 010 and 011
  node010->trainPrediction_ = 3.99;
  node011->trainPrediction_ = 4.3;

  // Children of node11
  Node::Node* node110 = node11->leftChild_;
  Node::Node* node111 = node11->rightChild_;

  // Set nodes 110 and 111
  node110->trainPrediction_ = 6.5;
  node111->trainPrediction_ = 7.1;


  // OK so now that we've constructed a small decision tree, let's put it to use! 

  

  delete partitionSequence;

}

// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( RootNodeTest ); 

#endif // ROOTNODETEST_HPP
