#include <cctype>
#include <iostream>
#include <cassert>
#include <sstream>
#include "argparse.hpp"

ArgParse::ArgParse(const int argc, char* argv[])
{

  tokens_.clear();

  if(argc == 1)
    {
      return;
    }

  int iter = 1;
  //pair<set<Token>::iterator, bool> lastAdded;
  while(iter < argc)
    {
      //First read the key(s)
      string str = argv[iter];
      vector<string> keys;
      if(getKeys(str,keys))
	{
	  assert(keys.size() > 0);
	  //cout << "keys: ";
	  for(size_t j = 0; j < keys.size(); ++j)
	    {
	      //cout << " '" << keys[j] << "'";
	      Token token;
	      token.name = keys[j];
	      token.value = "";
	      tokens_.push_back(token);
	    }
	  //cout << endl;
	}
      else
	{
	  tokens_[ tokens_.size()-1 ].value = str;
	}
      ++iter;
    }

  //printTokens();

  assert(!foundRepeats());

}

ArgParse::~ArgParse()
{

}

bool ArgParse::foundRepeats()
{
  string nameCopy;
  for(size_t i = 0; i < tokens_.size(); ++i)
    {
      bool foundOne = false;
      nameCopy = tokens_[i].name;
      for(size_t j = i; j < tokens_.size(); ++j)
	{
	  if(tokens_[j].name.compare(nameCopy) == 0)
	    {
	      if(!foundOne)
		{
		  foundOne = true;
		}
	      else
		{
		  return(true);
		}
	    }
	}
    }
  return(false);
}

/*
  template <typename T> void ArgParse::getArgument(const string& shortName, const string& longName, T& argument)
  {
  bool foundShort = false;
  bool foundLong = false;
  size_t foundAtIdx;
  for(size_t i = 0; i < tokens_.size(); ++i)
  {
  
  if(tokens_[i].name.compare(shortName) == 0)
  {
  foundShort = true;
  foundAtIdx = i;
  }
  
  if(tokens_[i].name.compare(longName) == 0)
  {
  foundLong = true;
  foundAtIdx = i;
  }
  }
  
  if(!foundShort && !foundLong)
  {
  argument = false;
  return;
  }
  
  assert(foundShort^foundLong);
  
  if(tokens_[foundAtIdx].value == "")
  {
  argument = true;
  return;
  }
  
  stringstream ss(tokens_[foundAtIdx].value);
  ss >> argument;
  assert(ss.eof());
  
  }
*/

bool ArgParse::getKeys(string& str, vector<string>& keys)
{
  keys.clear();

  //cout << "interpreting '" << str << "'..." << endl;

  //If the key starts as a short...
  if(str.size() > 1 && str.compare(0,1,"-") == 0 && str.compare(1,1,"-") != 0)
    {
      //... it's also possible that it's a bundle of shorts...
      for(size_t i = 1; i < str.size(); ++i)
	{
	  if(!isalpha(str[i]))
	    {
	      keys.clear();
	      return(false);
	    }
	  string newKey(str.substr(i,1));
	  keys.push_back(newKey);
	}

      return(true);
    }
  //If the key starts as a long...
  else if(str.size() > 2 && str.compare(0,2,"--") == 0)
    {
      for(size_t i = 2; i < str.size(); ++i)
	{
	  if(!isalpha(str[i]))
	    {
	      keys.clear();
	      return(false);
	    }
	}
      string newKey(str.substr(2));
      keys.push_back(newKey);
      return(true);
    }
  
  return(false);

}

void ArgParse::printTokens()
{
  for(size_t i = 0; i < tokens_.size(); ++i)
    {
      cout << tokens_[i].name << ": ";
      if(tokens_[i].value == "")
	{
	  cout << "true" << endl;
	}
      else
	{
	  cout << tokens_[i].value << endl;
	}
    }
}
