#include "utils.hpp"

#include <stdio.h>
#include <string.h>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <string>
#include <ios>

#include "murmurhash3.hpp"
#include "math.hpp"

string utils::tolower(const string& str) {

  string strcopy(str);
  transform(strcopy.begin(),strcopy.end(),strcopy.begin(),::tolower);
  return(strcopy);
}

// Returns a copy of input vector x with NAN-entries removed
// NOTE: is just a wrapper of the algorithm "remove_if"
vector<num_t> utils::removeNANs(vector<num_t> x) {

  x.erase( remove_if(x.begin(),x.end(),&datadefs::isNAN), x.end() );

  return( x );
}

string utils::num2str(const num_t x) {

  if ( datadefs::isNAN(x) ) return( datadefs::STR_NAN );
  
  stringstream ss;
  ss << x;
  
  return( ss.str() );
  
}

void utils::strv2numv(const vector<string>& strvec,
                         vector<datadefs::num_t>& numvec) {
  size_t n = strvec.size();
  numvec.resize(n);
  
  for(size_t strIdx = 0; strIdx < n; ++strIdx) {
    numvec[strIdx] = utils::str2<datadefs::num_t>(strvec[strIdx]);
  }
}

void utils::strv2catv(const vector<string>& strvec,
		      vector<datadefs::num_t>& catvec,
		      map<string,datadefs::num_t>& mapping,
		      map<datadefs::num_t,string>& backMapping) {
  
  size_t n = strvec.size();
  catvec.resize(n);

  mapping.clear();
  backMapping.clear();

  num_t val = static_cast<datadefs::num_t>(0);

  //Map unique strings to values and store values in catvec as doubles
  for(size_t strIdx = 0; strIdx < n; ++strIdx) {

    //If the string is not NaN ...
    if(!datadefs::isNAN_STR(strvec[strIdx])) {
      map<string,num_t>::iterator it;

      //Try to find the string in the map. If it's not found, extend the map...
      it = mapping.find(strvec[strIdx]);
      if(it == mapping.end()) {
        mapping.insert(pair<string,num_t>(strvec[strIdx],val));
        backMapping.insert(pair<num_t,string>(val,strvec[strIdx]));
        catvec[strIdx] = val;
        val = static_cast<num_t>(mapping.size());
      } else {
        catvec[strIdx] = it->second;
      }

    } else {    //If the string is defined to NaN, however...
      catvec[strIdx] = datadefs::NUM_NAN;
    }
  }

}

void utils::sortDataAndMakeRef(const bool isIncreasingOrder,
			       vector<num_t>& data,
			       vector<size_t>& refIcs) {

  //assert(v.size() == ref_ics.size());
  vector<pair<num_t,size_t> > pairedv(data.size()); // !! Understandibility:
                                                    // !! consider a typedef
                                                    // !! pairedv, leaving the
                                                    // !! actual variable name
                                                    // !! as something more
                                                    // !! descriptive.

  refIcs = utils::range(data.size());

  datadefs::make_pairedv<num_t,size_t>(data,refIcs,pairedv);

  //pairedv.erase(remove_if(pairedv.begin(),pairedv.end(),&datadefs::pairedIsNAN), pairedv.end());

  if(isIncreasingOrder) {
    sort(pairedv.begin(),pairedv.end(),datadefs::increasingOrder<size_t>());
  } else {
    sort(pairedv.begin(),pairedv.end(),datadefs::decreasingOrder<size_t>());
  }

  datadefs::separate_pairedv<num_t,size_t>(pairedv,data,refIcs);
}


// Removes all newline and any trailing characters
string utils::chomp(const string& str, const string& nl) {
  
  size_t endStr = str.find_first_of(nl);
  return( str.substr(0, endStr) );

}

// Remove all leading and trailing whitespaces
string utils::trim(const string& str, const string& wh) {
  
  size_t beginStr = str.find_first_not_of(wh);
  if ( beginStr == string::npos ) {
    // no content
    return("");
  }

  size_t endStr = str.find_last_not_of(wh);
  size_t range = endStr - beginStr + 1;
  
  return( str.substr(beginStr, range) );
}

set<string> utils::keys(const string& str, const char delimiter) {

  set<string> ret;

  vector<string> items = utils::split(str,delimiter);

  for ( size_t i = 0; i < items.size(); ++i ) {
    ret.insert(items[i]);
  }
  
  return( ret );

}

map<string,string> utils::parse(const string& str,
				const char delimiter,
				const char separator,
				const char comment) {

  stringstream streamObj(str);

  return( utils::parse(streamObj,delimiter,separator,comment) );

}

map<string,string> utils::parse(istream& streamObj,
				const char delimiter,
				const char separator,
				const char comment) {

  map<string,string> ret;

  string key;

  while ( !streamObj.eof() ) {

    // Parse the key
    getline(streamObj,key,separator);
    assert( streamObj.good() );
    assert( !streamObj.eof() );
    
    // Peek the next characeter and check if it's a comment 
    if ( streamObj.peek() == comment ) {

      // ignore the comment character...
      streamObj.ignore();

      // ... and get the value for the key
      getline(streamObj,ret[key],comment);

      assert( ret.find(key) != ret.end() );

      // If the next character is a delimiter, ignore it
      if ( streamObj.peek() == delimiter ) {
	streamObj.ignore();
      } 

    } else {
      
      // Without the comment character we just read until the next delimiter
      getline(streamObj,ret[key],delimiter);
      
    }
 
  }
  
  // The possible carriage return and end-of-line characters need to be removed
  ret[key] = utils::chomp(ret[key]);
    
  return( ret );

}

vector<string> utils::split(const string& str, const char delimiter, const string& wh) {
  stringstream streamObj(str);
  return( utils::split(streamObj,delimiter,wh) );
}

vector<string> utils::split(istream& streamObj, const char delimiter, const string& wh) {

  string newItem("");
  vector<string> items;

  while ( getline(streamObj,newItem,delimiter) ) {
    newItem = utils::trim(newItem,wh);
    items.push_back(newItem);
  }

  return( items );

}

unordered_set<uint32_t> utils::hashText(const string& text) {

  unordered_set<uint32_t> hashes;

  char const* p = text.c_str();
  char const* q = strpbrk(p+1,datadefs::tokenDelimiters);
  for ( ; ; q = strpbrk(p,datadefs::tokenDelimiters) ) {
    if ( distance(p,q) > 0 ) {
      string token(p,q);
      uint32_t h;
      MurmurHash3_x86_32(utils::tolower(token).c_str(),token.length(),0,&h);
      hashes.insert( h );
    }
    if ( q != NULL ) {
      p = q + 1;
    } else {
      string token2(p);
      uint32_t h;
      MurmurHash3_x86_32(utils::tolower(token2).c_str(),token2.length(),0,&h);
      hashes.insert( h );
      break;
    }
  }
   
  return(hashes);
}

vector<string> utils::readListFromFile(const string& fileName, const char delimiter) {
  
  ifstream streamObj( fileName.c_str() );

  if ( ! streamObj.good() ) {
    cerr << "ERROR: file '" << fileName << "' could not be opened for reading. Check that the file exists and is not corrupted." << endl;
    exit(1);
  }
  
  return( utils::split(streamObj,delimiter) );
}

void utils::filterSort(const bool isIncreasingOrder,
		       vector<num_t>& data,
		       vector<size_t>& refIcs) {
  
  //assert(v.size() == ref_ics.size());
  vector<pair<num_t,size_t> > pairedv(data.size()); // !! Understandibility:
  // !! consider a typedef
  // !! pairedv, leaving the
  // !! actual variable name
  // !! as something more
  // !! descriptive.
  refIcs = utils::range( data.size() );
  
  datadefs::make_pairedv<num_t,size_t>(data,refIcs,pairedv);
  
  pairedv.erase(remove_if(pairedv.begin(),pairedv.end(),&datadefs::pairedIsNAN), pairedv.end());
  
  if(isIncreasingOrder) {
    sort(pairedv.begin(),pairedv.end(),datadefs::increasingOrder<size_t>());
  } else {
    sort(pairedv.begin(),pairedv.end(),datadefs::decreasingOrder<size_t>());
  }
  
  datadefs::separate_pairedv<num_t,size_t>(pairedv,data,refIcs);
}

// !! Documentation: this is just a drop-in replacement for Python's range()
// !! function, hinted by the size of the input vector. It mutates ics,
// !! specifying a 0-based range in all cases, and could be made more robust if
// !! the starting value could be given.
vector<size_t> utils::range(const size_t n) {
  
  vector<size_t> ics(n);
  
  for(size_t i = 0; i < n; ++i) {
    ics[i] = i;
  }

  return( ics );

}


istream& utils::safeGetline(istream& is, string& t) {

  t.clear();

  // The characters in the stream are read one-by-one using a std::streambuf.
  // That is faster than reading them one-by-one using the std::istream.
  // Code that uses streambuf this way must be guarded by a sentry object.
  // The sentry object performs various tasks,
  // such as thread synchronization and updating the stream state.

  istream::sentry se(is, true);
  streambuf* sb = is.rdbuf();

  for(;;) {
    int c = sb->sbumpc();
    switch (c) {
    case '\r':
      c = sb->sgetc();
      if(c == '\n')
	sb->sbumpc();
      return is;
    case '\n':
    case EOF:
      return is;
    default:
      t += (char)c;
    }
  }
}

vector<vector<size_t> > utils::splitRange(const size_t nElements, const size_t nSplits) {

  vector<size_t> ics = utils::range(nElements);

  assert( nSplits >= 1 );

  vector<vector<size_t> > splits(nSplits);

  size_t splitSize = ics.size() / nSplits ;

  size_t remainder = ics.size() % nSplits;

  size_t startIdx = 0;

  for ( size_t splitIdx = 0; splitIdx < nSplits; ++splitIdx ) {

    size_t stopIdx = startIdx + splitSize;

    if ( splitIdx < remainder ) stopIdx++;

    splits[splitIdx].resize(stopIdx-startIdx);

    copy(ics.begin()+startIdx,ics.begin()+stopIdx,splits[splitIdx].begin());

    startIdx = stopIdx;

  }

  return( splits );

}

num_t utils::numericalFeatureSplitsNumericalTarget(const vector<num_t>& tv,
						   const vector<num_t>& fv,
						   const size_t minSamples,
						   size_t& splitIdx) {

  size_t n_tot = tv.size();
  size_t n_left = 0;
  size_t n_right = n_tot;

  // We start with all samples on the left branch
  num_t mu_tot = math::mean(tv);
  num_t mu_left = 0.0;
  num_t mu_right = mu_tot;

  // Make sure the squared error didn't become corrupted by NANs
  assert( !datadefs::isNAN(mu_tot) );

  num_t DI_best = 0.0;

  // Add samples one by one from left to right until we hit the
  // minimum allowed size of the branch
  for( size_t i = 0; i < n_tot - minSamples; ++i ) {

    // Add n'th sample tv[i] from left to right and update
    // mean and squared error
    ++n_left;
    --n_right;
    mu_left  += ( tv[i] - mu_left  ) / n_left;
    mu_right -= ( tv[i] - mu_right ) / n_right;

    // If the sample is repeated and we can continue, continue
    if ( n_left < minSamples || (n_left < n_tot - minSamples && fv[ i + 1 ] == fv[ i ]) ) {
      continue;
    }

    // If the current split point yields a better split than the best,
    // update DI_best and bestSplitIdx
    num_t DI = math::deltaImpurity_regr(mu_tot,n_tot,mu_left,n_left,mu_right,n_right);
    //cout << "(" << n_left << "," << tv[i] << "," << DI << ")";
    if (  DI > DI_best ) {

      splitIdx = i;
      DI_best = DI;

    }

  }

  return(DI_best);
}

num_t utils::numericalFeatureSplitsCategoricalTarget(const vector<num_t>& tv,
						     const vector<num_t>& fv,
						     const size_t minSamples,
						     size_t& splitIdx) {
 
  size_t n_tot = tv.size();
  size_t n_left = 0;
  size_t n_right = n_tot;
 
  map<num_t,size_t> freq_right;
  size_t sf_right = 0;
  
  for ( size_t i = 0; i < n_tot; ++i ) {
    math::incrementSquaredFrequency(tv[i],freq_right,sf_right);
  }
  
  size_t sf_tot = sf_right;
  
  map<num_t,size_t> freq_left;
  size_t sf_left = 0;

  num_t DI_best = 0.0;
  
  // Add samples one by one from right to left until we hit the
  // minimum allowed size of the branch
  for( size_t i = 0; i < n_tot - minSamples; ++i ) {
    
    // Add n'th sample tv[i] from right to left and update
    // mean and squared frequency
    math::incrementSquaredFrequency(tv[i],freq_left,sf_left);
    ++n_left;
    
    math::decrementSquaredFrequency(tv[i],freq_right,sf_right);
    --n_right;
    
    // If we have repeated samples and can continue, continue
    if ( n_left < minSamples || (n_left < n_tot - minSamples && fv[ i + 1 ] == fv[ i ]) ) {
      continue;
    }
    
    // If the split point "i-1" yields a better split than the previous one,
    // update se_best and bestSplitIdx
    num_t DI = math::deltaImpurity_class(sf_tot,n_tot,sf_left,n_left,sf_right,n_right);
    //cout << " tv=" << features_[targetIdx].backMapping[tv[i]] << " fv=" << fv[i] << " nl=" << n_left << " sfl=" << sf_left << " nr=" << n_right << " sfr=" << sf_right << " nt=" << n_tot << " sft=" << sf_tot << " DI=" << DI << endl;
    
    if ( DI > DI_best ) {
      splitIdx = i;
      DI_best = DI;
    }
    
  }
  
  return(DI_best);
  
}

num_t utils::categoricalFeatureSplitsNumericalTarget(const vector<num_t>& tv,
						     const vector<num_t>& fv,
						     const size_t minSamples,
						     map<num_t,vector<size_t> >& fmap_left,
						     map<num_t,vector<size_t> >& fmap_right) {
  

  fmap_left.clear();

  size_t n_tot = 0;
  datadefs::map_data(fv,fmap_right,n_tot);
  size_t n_right = n_tot;
  size_t n_left = 0;

  num_t mu_tot = math::mean(tv);
  num_t mu_right = mu_tot;
  num_t mu_left = 0.0;

  num_t DI_best = 0.0;

  while ( fmap_right.size() > 1 ) {

    map<num_t,vector<size_t> >::iterator it_best( fmap_right.end() );

    // We test each category one by one and see if the fitness becomes improved
    for ( map<num_t,vector<size_t> >::iterator it( fmap_right.begin() ); it != fmap_right.end() ; ++it ) {

      //cout << "Testing to split with feature '" << treedata->getRawFeatureData(featureIdx,it->first) << "'" << endl;

      // Take samples from right and put them left
      //cout << "from right to left: [";
      size_t n_left_c = n_left;
      size_t n_right_c = n_right;
      num_t mu_left_c = mu_left;
      num_t mu_right_c = mu_right;

      for(size_t i = 0; i < it->second.size(); ++i) {
	//cout << " " << it->second[i];

	++n_left_c;
	--n_right_c;
	mu_left_c  += ( tv[ it->second[i] ] - mu_left_c  ) / n_left_c;
	mu_right_c -= ( tv[ it->second[i] ] - mu_right_c ) / n_right_c;

      }
      //cout << " ]" << endl;

      //If the split reduces impurity even further, save the point
      num_t DI = math::deltaImpurity_regr(mu_tot,n_tot,mu_left_c,n_left_c,mu_right_c,n_right_c);
      if ( DI > DI_best ) { //&& n_left >= minSamples && n_right >= minSamples )

	it_best = it;
	DI_best = DI;
      }

    }

    // After testing all categories,
    // if we couldn't find any split that would reduce impurity,
    // we'll exit the loop
    if ( it_best == fmap_right.end() ) {
      //cout << " -- STOP --" << endl;
      break;
    }

    // Otherwise move samples from right to left
    for(size_t i = 0; i < it_best->second.size(); ++i) {
      //cout << " " << it->second[i];

      ++n_left;
      --n_right;
      mu_left  += ( tv[ it_best->second[i] ] - mu_left  ) / n_left;
      mu_right -= ( tv[ it_best->second[i] ] - mu_right ) / n_right;

    }

    // Update the maps
    fmap_left.insert( *it_best );
    fmap_right.erase( it_best->first );

  }

  assert( n_left + n_right == n_tot );

  if( n_left < minSamples || n_right < minSamples ) {
    DI_best = 0.0;
  }

  return(DI_best);

}

num_t utils::categoricalFeatureSplitsCategoricalTarget(const vector<num_t>& tv,
						       const vector<num_t>& fv,
						       const size_t minSamples,
						       map<num_t,vector<size_t> >& fmap_left,
						       map<num_t,vector<size_t> >& fmap_right) {
  
  fmap_left.clear();

  size_t n_tot = 0;
  datadefs::map_data(fv,fmap_right,n_tot);
  size_t n_right = n_tot;
  size_t n_left = 0;

  map<num_t,size_t> freq_left,freq_right;
  size_t sf_left = 0;
  size_t sf_right = 0;

  for( size_t i = 0; i < n_tot; ++i ) {
    math::incrementSquaredFrequency(tv[i], freq_right, sf_right);
  }

  size_t sf_tot = sf_right;

  num_t DI_best = 0.0;

  while ( fmap_right.size() > 1 ) {

    map<num_t,vector<size_t> >::iterator it_best( fmap_right.end() );
    //cout << "There are " << fmap_right.size() << " categories on right" << endl;

    // We test each category one by one and see if the fitness becomes improved
    for ( map<num_t,vector<size_t> >::iterator it( fmap_right.begin() ); it != fmap_right.end() ; ++it ) {

      //cout << "Testing to split with feature '" << treedata->getRawFeatureData(featureIdx,it->first) << "'" << endl;

      // Take samples from right and put them left
      //cout << "from right to left: [";
      for(size_t i = 0; i < it->second.size(); ++i) {

	// Add sample to left
	++n_left;
	math::incrementSquaredFrequency(tv[ it->second[i] ], freq_left, sf_left);

	// Remove sample from right
	--n_right;
	math::decrementSquaredFrequency(tv[ it->second[i] ], freq_right, sf_right);

      }
      //cout << " ]" << endl;

      //If the impurity becomes reduced even further, save the point
      num_t DI = math::deltaImpurity_class(sf_tot,n_tot,sf_left,n_left,sf_right,n_right);
      if ( DI > DI_best ) { //&& n_left >= minSamples && n_right >= minSamples )

	it_best = it;
	DI_best = DI;
      }

      // Take samples from left and put them right
      //cout << "From left to right: [";
      for(size_t i = 0; i < it->second.size(); ++i) {

	// Add sample to right
	++n_right;
	math::incrementSquaredFrequency(tv[ it->second[i] ], freq_right, sf_right);

	// Remove sample from left
	--n_left;
	math::decrementSquaredFrequency(tv[ it->second[i] ], freq_left, sf_left);

      }
      //cout << " ]" << endl;

    }

    // After testing all categories,
    // if we couldn't find any split that would reduce impurity,
    // we'll exit the loop
    if ( it_best == fmap_right.end() ) {
      //cout << " -- STOP --" << endl;
      break;
    }

    // Take samples from right and put them left
    for(size_t i = 0; i < it_best->second.size(); ++i) {

      // Add sample to left
      ++n_left;
      math::incrementSquaredFrequency(tv[ it_best->second[i] ], freq_left, sf_left);

      // Remove sample from right
      --n_right;
      math::decrementSquaredFrequency(tv[ it_best->second[i] ], freq_right, sf_right);

    }

    // Update the maps
    fmap_left.insert( *it_best );
    fmap_right.erase( it_best->first );

  }

  assert( n_left + n_right == n_tot );

  if( n_left < minSamples || n_right < minSamples ) {
    DI_best = 0.0;
  }

  return(DI_best);

}
