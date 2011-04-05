#ifndef GBT_HPP
#define GBT_HPP

#include <cstdlib>
#include "node.hpp"
#include "treedata.hpp"

using namespace std;

class GBT
{
public:

	// Constructs Gradient Boosting Tree "forest".
	// nmaxleaves is typically a small value, such as 6.
	// shrinkage is between 0 and 1, typically something like 0.125
	// subSampleSize is typically 0.4 to 0.6
	GBT(Treedata* treeData, size_t nTrees, size_t nMaxLeaves, float shrinkage, float subSampleSize);
	~GBT();

	// actual growing
	void growForest();
	// produce predictions of a data set
	void predictForest(Treedata* treeData, vector<num_t>& prediction);

	//Selects the target feature that is to be predicted
	void selectTarget(size_t targetidx);
	//Gets the selected target feature
	size_t getTarget();

private:

	void growTree(size_t treeidx);
	num_t predictSampleByTree(size_t sampleidx, size_t treeidx, Treedata* treeData);
	void predictDatasetByTree(size_t treeidx, vector<num_t>& prediction);
	void recursiveNodesplit(size_t treeidx, size_t nodeidx, vector<size_t>& sampleics);
	void SetNodePrediction( size_t treeidx, size_t nodeidx, vector<size_t>& sampleics);

	Treedata* treeData_;
	size_t nTrees_;
	size_t nMaxLeaves_;
	size_t nMaxNodes_;
	float shrinkage_;
	float subSampleSize_;

	vector<size_t> nNodes_; //Number of used nodes in each tree.
	vector<vector<Node> > forest_; // forest_[i][j] is the j'th node of i'th tree. forest_[i][0] is the rootnode.
};

#endif
