#include <cstdlib>
#include <cassert>
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include "randomforest.hpp"
//#include "GBT.hpp"
#include "treedata.hpp"
#include "datadefs.hpp"

using namespace std;
using datadefs::num_t;

namespace po = boost::program_options;

enum MatrixFormat { FEATURE_ROWS, FEATURE_COLUMNS };

const size_t DEFAULT_NODESIZE = 5;
const size_t DEFAULT_NPERMS = 9;
const num_t DEFAULT_PTHRESHOLD = 0.5;

int main(int argc, char* argv[])
{
    try {
        // Declare the supported options.
        po::options_description general_opts("General options");
        general_opts.add_options()
          ("input,I", po::value<std::string>(), "Input feature matrix")
          ("format,f", po::value<std::string>(), "Feature matrix file format")
          ("targetidx,i", po::value<size_t>(), "Index of target feature (0-base)")
	  ("output,O", po::value<std::string>(), "Output association file")
        ;
        
        po::options_description rf_opts("Random forest options");
        rf_opts.add_options()
          ("ntrees,n", po::value<size_t>(), "Number of decision trees")
          ("mtry,m", po::value<size_t>(), "Number of randomly selected features for each split process")
          ("nodesize,s", po::value<size_t>(), "Minimum number of samples per node")
	  ("nperms,p", po::value<size_t>(), "Number of permutations")
	  ("pthreshold,t", po::value<num_t>(), "p-value threshold for output associations")
          ("help,h", "Display help message")
        ;
        
        po::options_description all_options("Allowed options");
        all_options.add(general_opts).add(rf_opts);
        
        po::variables_map var_map;
        po::store(po::parse_command_line(argc, argv, all_options), var_map);
        po::notify(var_map);    
        
        // Print help if needed and quit.
        if (var_map.count("help")) {
            cout << all_options << "\n";
            return(EXIT_FAILURE);
        }
        
        // Parse feature matrix file name
        std::string matrix_filename;
        if (var_map.count("input")) {
            matrix_filename = var_map["input"].as<std::string>();
            cout << "Matrix: " << matrix_filename << endl;
        } else {
            cout << "Feature matrix file name not set, quitting." << endl;
            cout << all_options << "\n";
            return(EXIT_FAILURE);
        }    

        // Parse output file name
	std::string output_filename;
        if (var_map.count("output")) {
	  output_filename = var_map["output"].as<std::string>();
	  cout << "Output: " << output_filename << endl;
        } else {
	  cout << "Output file name not set, quitting." << endl;
	  cout << all_options << "\n";
	  return(EXIT_FAILURE);
        }
        
        // Parse feature matrix file format
        enum MatrixFormat mat_format = FEATURE_ROWS;
        
        if (var_map.count("format")) {
            std::string format_string = var_map["format"].as<std::string>();
            if (boost::iequals(format_string, "frows")) {
                mat_format = FEATURE_ROWS;
            } else if (boost::iequals(format_string, "fcols")) {
                mat_format = FEATURE_COLUMNS;
            } else {
                cout << "Unknown feature matrix file format \"" << format_string << "\", quitting." << endl;
                cout << all_options << "\n";
                return(EXIT_FAILURE);
            }
        } else {
            cout << "Feature matrix file format not set, defaulting to FEATURE_ROWS." << endl;
            mat_format = FEATURE_ROWS;
        }
        
        size_t targetidx = 0;
        if (var_map.count("targetidx")) {
            targetidx = var_map["targetidx"].as<size_t>();
        } else {
            cout << "Parameter \"targetidx\" not set, quitting." << endl;
            cout << all_options << "\n";
            return(EXIT_FAILURE);
        }
        
        size_t ntrees = 0;
        if (var_map.count("ntrees")) {
            ntrees = var_map["ntrees"].as<size_t>();
        } else {
            cout << "Parameter 'ntrees' not set, defaulting to number of samples in the feature matrix." << endl;
        }
        
        size_t mtry = 0;
        if (var_map.count("mtry")) {
            mtry = var_map["mtry"].as<size_t>();
        } else {
            cout << "Parameter 'mtry' not set, defaulting to number of 10% of features in the feature matrix." << endl;
        }
        
        size_t nodesize = 0;
        if (var_map.count("nodesize")) {
            nodesize = var_map["nodesize"].as<size_t>();
        } else {
            cout << "Parameter 'nodesize' not set, defaulting to " << DEFAULT_NODESIZE << "." << endl;
	    nodesize = DEFAULT_NODESIZE;
        }
        
        size_t nperms = 0;
        if (var_map.count("nperms")) {
	  nperms = var_map["nperms"].as<size_t>();
        } else {
	  cout << "Parameter 'nperms' not set, defaulting to " << DEFAULT_NPERMS << "." << endl;
	  nperms = DEFAULT_NPERMS;
        }

        num_t pthreshold = 0;
        if (var_map.count("pthreshold")) {
          pthreshold = var_map["pthreshold"].as<num_t>();
        } else {
          cout << "Parameter 'pthreshold' not set, defaulting to " << DEFAULT_PTHRESHOLD << "." << endl;
          pthreshold = DEFAULT_PTHRESHOLD;
        }



        bool is_featurerows = true;
        if (mat_format == FEATURE_COLUMNS) {
            is_featurerows = false;
        }
        
	//Read data into Treedata object
        Treedata treedata(matrix_filename, is_featurerows);
        
	assert(treedata.nfeatures() >= mtry);
	assert(treedata.nsamples() > 2*nodesize);

        //Construct a Random Forest object
        Randomforest RF(&treedata,ntrees,mtry,nodesize);
	RF.select_target(targetidx);
	//treedata.print();   
	
	//size_t nperms = 9;
	num_t alpha = 0.50;
        vector<num_t> pvalues(treedata.nfeatures());

	//clock_t time_start(clock());
	RF.grow_forest(nperms,alpha,pvalues);
	//cout << "Time elapsed: " << float(clock() - time_start)/CLOCKS_PER_SEC << " seconds" << endl;

	vector<size_t> ref_ics(treedata.nfeatures());
	//vector<string> fnames = treedata.featureheaders();
	datadefs::sort_and_make_ref<num_t>(pvalues,ref_ics);
	//datadefs::sort_from_ref<string>(fnames,ref_ics);
	
	string target_str = treedata.get_featureheader(targetidx);

	ofstream os(output_filename.c_str());
	for(size_t i = 0; i < treedata.nfeatures(); ++i)
	  {
	    if(pvalues[i] > pthreshold)
	      {
		break;
	      }
	    os << target_str << "\t" << treedata.get_featureheader(ref_ics[i]) << "\t" 
	       << pvalues[i] << "\t" << treedata.corr(targetidx,ref_ics[i]) << endl;
	  }
	os.close();
	
    }
    catch(exception& e) {
        cerr << e.what() << "\n";
    }

    return(EXIT_SUCCESS);
}
