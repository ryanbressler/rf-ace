#include "treedata.hpp"
#include<cstdlib>
#include<fstream>

using namespace std;

Treedata::Treedata(string file_str, bool is_featurerows)
{
  if(is_featurerows)
    {
      this->read_featurerows();
    }
  else
    {
      this->read_featurecols();
    }
}

Treedata::~Treedata()
{
}

void Treedata::read_featurerows()
{

}

void Treedata::read_featurecols()
{

}


