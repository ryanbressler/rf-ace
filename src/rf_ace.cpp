#include <cstdlib>
#include <cassert>
#include <string>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include "randomforest.hpp"
#include "treedata.hpp"

using namespace std;

namespace po = boost::program_options;

enum MatrixFormat { FEATURE_ROWS, FEATURE_COLUMNS };

const size_t DEFAULT_NODESIZE = 5;

int main(int argc, char* argv[])
{
    try {
        // Declare the supported options.
        po::options_description general_opts("General options");
        general_opts.add_options()
          ("filename", po::value<std::string>(), "Name of feature matrix file")
          ("format", po::value<std::string>(), "Featrure matrix file format")
          ("targetidx,i", po::value<size_t>(), "Index of target feature")
        ;
        
        po::options_description rf_opts("Random forest options");
        rf_opts.add_options()
          ("ntrees,n", po::value<size_t>(), "Number of decision trees")
          ("mtry,m", po::value<size_t>(), "Number of randomly selected features for each split process")
          ("nodesize,s", po::value<size_t>(), "Minimum number of samples per node")
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
        if (var_map.count("filename")) {
            matrix_filename = var_map["filename"].as<std::string>();
            cout << "Matrix: " << matrix_filename << endl;
        } else {
            cout << "Feature matrix file name not set, quitting." << endl;
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
        }
        
        //FIRST PART: read data into Treedata class (features are rows)
        bool is_featurerows = true;
        if (mat_format == FEATURE_COLUMNS) {
            is_featurerows = false;
        }
        
        Treedata treedata(matrix_filename, is_featurerows);
        
	assert(treedata.nfeatures() >= mtry);
	assert(treedata.nsamples() > 2*nodesize);

        //SECOND PART: construct a Random Forest object
        Randomforest RF(&treedata,ntrees,mtry,nodesize);
        
	RF.select_target(targetidx);
	treedata.print();
	
	RF.grow_forest();
        
        return(EXIT_SUCCESS);
    }
    catch(exception& e) {
        cerr << e.what() << "\n";
    }    
}
