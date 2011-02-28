//treedata.hpp
//
//

#include<cstdlib>

using namespace std;

class Treedata 
{
public:
  Treedata(string file_str,bool is_feature_rows);
  ~Treedata();
  
private:

  void read_featurerows();
  void read_featurecols();

};
