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
  //CPPUNIT_TEST( test_getSplitter );
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
  //CPPUNIT_TEST( test_splitFitness );
  CPPUNIT_TEST( test_recursiveNDescendantNodes );
  CPPUNIT_TEST_SUITE_END();
  
public:
  void setUp();
  void tearDown();

  void test_setSplitter();
  //void test_getSplitter();
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
  //void test_splitFitness();
  void test_recursiveNDescendantNodes();

};

void NodeTest::setUp() {}
void NodeTest::tearDown() {}


void NodeTest::test_setSplitter() {
  
  size_t splitterIdx = 3;
  datadefs::num_t splitLeftLeqValue = 0.5; 
  
  //Splitter::Splitter splitter(0.5);
  
  Node::Node node;
  
  //CPPUNIT_ASSERT( false );
  
  CPPUNIT_ASSERT( !node.splitter_ );
  
  node.setSplitter(splitterIdx,"foo",splitLeftLeqValue);
  
  CPPUNIT_ASSERT( node.splitterIdx_ == splitterIdx );
  CPPUNIT_ASSERT( node.splitter_->splitterType_ == Splitter::NUMERICAL_SPLITTER );
  CPPUNIT_ASSERT( fabs(node.splitter_->splitLeftLeqValue_ - splitLeftLeqValue) < datadefs::EPS );
  
  
}


//void NodeTest::test_getSplitter() {
  //Node::getSplitter(...);
//}

void NodeTest::test_percolateData() {
  
  Node::Node node0;
  //Splitter splitter("foo",0.1);
  node0.setSplitter(1,"foo",0.1);
  CPPUNIT_ASSERT( node0.leftChild_ == node0.percolateData(0.09) );
  CPPUNIT_ASSERT( node0.rightChild_ == node0.percolateData(0.11) );

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

  Treedata treedata("test_2by8_numerical_matrix.tsv",'\t',':');

  CPPUNIT_ASSERT( treedata.getFeatureName(0) == "N:F1" );
  CPPUNIT_ASSERT( treedata.getFeatureName(1) == "N:F2" );

  vector<size_t> sampleIcs_left(0);
  vector<size_t> sampleIcs_right(8);
  datadefs::range(sampleIcs_right);
  
  datadefs::num_t splitValue;
  datadefs::num_t splitFitness;

  Node node;

  Node::GrowInstructions GI;
  GI.minNodeSizeToStop = 2;

  size_t targetIdx = 1;
  size_t featureIdx = 0;

  node.numericalFeatureSplit(&treedata,
			     targetIdx,
			     featureIdx,
			     GI,
			     sampleIcs_left,
			     sampleIcs_right,
			     splitValue,
			     splitFitness);

  //cout << splitValue << endl;

  Splitter::Splitter splitter("foo",splitValue);

  CPPUNIT_ASSERT( splitter.splitsLeft(3.1 - 0.00001) );
  CPPUNIT_ASSERT( splitter.splitsRight(3.1) );

  CPPUNIT_ASSERT( fabs(splitFitness - 0.530492285084497) < 1e-10 );

  //datadefs::print(sampleIcs_left);
  //datadefs::print(sampleIcs_right);

  CPPUNIT_ASSERT(sampleIcs_left.size() == 4);
  for(size_t i = 0; i < sampleIcs_left.size(); ++i) {
    size_t idx = sampleIcs_left[i];
    CPPUNIT_ASSERT(idx == 1 || idx == 6 || idx == 0 || idx == 7);
  }

  CPPUNIT_ASSERT(sampleIcs_right.size() == 2);
  for(size_t i = 0; i < sampleIcs_right.size(); ++i) {
    size_t idx = sampleIcs_right[i];
    CPPUNIT_ASSERT(idx == 5 || idx == 4);
  }

}

void NodeTest::test_categoricalFeatureSplit() { 

  size_t n = 10;

  Treedata treedata("test_3by10_categorical_matrix.tsv",'\t',':');

  size_t featureIdx = 0;

  size_t targetIdx1 = 1; 
  size_t targetIdx2 = 2;

  Node node;
  PartitionSequence PS;

  Node::GrowInstructions GI;
  GI.minNodeSizeToStop = 2;
  GI.partitionSequence = &PS;

  vector<size_t> sampleIcs_left(0);
  vector<size_t> sampleIcs_right(n);
  datadefs::range(sampleIcs_right);
  
  set<num_t> splitValues_left,splitValues_right;
  datadefs::num_t splitFitness;
  
  node.categoricalFeatureSplit(&treedata,
			       targetIdx1,
			       featureIdx,
			       GI,
			       sampleIcs_left,
			       sampleIcs_right,
			       splitValues_left,
			       splitValues_right,
			       splitFitness);
   
  //datadefs::print<num_t>(splitValues_left);
  //datadefs::print<num_t>(splitValues_right);
  //datadefs::print<size_t>(sampleIcs_left);
  //datadefs::print<size_t>(sampleIcs_right);
  //datadefs::print<num_t>(treedata.getFeatureData(targetIdx1));
  //datadefs::print<num_t>(treedata.getFeatureData(featureIdx));
  //datadefs::print<num_t>(treedata.getFeatureData(targetIdx1,sampleIcs_left));
  //datadefs::print<num_t>(treedata.getFeatureData(targetIdx1,sampleIcs_right));
  
  set<string> rawSplitValues_left,rawSplitValues_right;

  for(set<num_t>::const_iterator it(splitValues_left.begin()); it != splitValues_left.end(); ++it ) {
    rawSplitValues_left.insert(treedata.getRawFeatureData(featureIdx,*it));
  }
  
  for(set<num_t>::const_iterator it(splitValues_right.begin()); it != splitValues_right.end(); ++it ) {
    rawSplitValues_right.insert(treedata.getRawFeatureData(featureIdx,*it));
  }

  Splitter::Splitter splitter("foo",rawSplitValues_left,rawSplitValues_right);
  
  CPPUNIT_ASSERT( sampleIcs_left.size() == sampleIcs_right.size() );
  
  CPPUNIT_ASSERT( splitter.splitsLeft("0") );
  CPPUNIT_ASSERT( splitter.splitsLeft("1") );
  CPPUNIT_ASSERT( splitter.splitsLeft("2") );
  CPPUNIT_ASSERT( splitter.splitsLeft("3") );
  CPPUNIT_ASSERT( splitter.splitsLeft("4") );

  //cout << iter << endl;
  
  for(size_t i = 0; i < sampleIcs_left.size(); ++i ) {
    //cout << " " << treedata.getRawFeatureData(targetIdx1,sampleIcs_left[i]); 
    CPPUNIT_ASSERT( treedata.getRawFeatureData(targetIdx1,sampleIcs_left[i]) == "1" );
    CPPUNIT_ASSERT( treedata.getRawFeatureData(targetIdx1,sampleIcs_right[i]) == "0" );
  }

  //cout << endl;
  
  //cout << splitFitness << endl;
  
  CPPUNIT_ASSERT( fabs(splitFitness - 1) < datadefs::EPS );
   
  sampleIcs_left.clear();
  sampleIcs_right.resize(n);
  datadefs::range(sampleIcs_right);
  
  splitValues_left.clear();
  splitValues_right.clear();
  //datadefs::num_t splitFitness;

  node.categoricalFeatureSplit(&treedata,
			       targetIdx2,
			       featureIdx,
			       GI,
			       sampleIcs_left,
			       sampleIcs_right,
			       splitValues_left,
			       splitValues_right,
			       splitFitness);
  
  //Splitter::Splitter splitter;
  //CPPUNIT_ASSERT( splitter.splitterType_ == Splitter::NO_SPLITTER );

  rawSplitValues_left.clear();
  rawSplitValues_right.clear();

  for(set<num_t>::const_iterator it(splitValues_left.begin()); it != splitValues_left.end(); ++it ) {
    rawSplitValues_left.insert(treedata.getRawFeatureData(featureIdx,*it));
  }

  for(set<num_t>::const_iterator it(splitValues_right.begin()); it != splitValues_right.end(); ++it ) {
    rawSplitValues_right.insert(treedata.getRawFeatureData(featureIdx,*it));
  }

  
  Splitter::Splitter splitter2("foo",rawSplitValues_left,rawSplitValues_right);
  
  CPPUNIT_ASSERT( splitter2.splitterType_ == Splitter::CATEGORICAL_SPLITTER );
  
  CPPUNIT_ASSERT( sampleIcs_left.size() == 4 );
  CPPUNIT_ASSERT( sampleIcs_right.size() == 6 );
  CPPUNIT_ASSERT( splitter2.splitsLeft("1") );
  CPPUNIT_ASSERT( splitter2.splitsLeft("2") );
  
  //if ( iter == 0 ) {
  CPPUNIT_ASSERT( fabs( splitFitness - 0.642857142857143 ) < 1e-10 );
  //} else {
  //CPPUNIT_ASSERT( fabs(splitFitness - 0.843750000000000 ) < 1e-10 );
  //}
  
}

//void NodeTest::test_splitFitness() { 
//}

void NodeTest::test_recursiveNDescendantNodes() {
  
}
// Registers the fixture into the test 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( NodeTest ); 

#endif // NODETEST_HPP
