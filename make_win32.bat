mkdir bin

cl /EHsc /O2 /analyze /Febin\rf-ace-filter-win32.exe src\rf_ace_filter.cpp src\statistics.cpp src\progress.cpp src\stochasticforest.cpp src\rootnode.cpp src\node.cpp src\splitter.cpp src\treedata.cpp src\mtrand.cpp src\datadefs.cpp src\math.cpp src\utils.cpp

cl /EHsc /O2 /analyze /Febin\rf-ace-build-predictor-win32.exe src\rf_ace_build_predictor.cpp src\statistics.cpp src\progress.cpp src\stochasticforest.cpp src\rootnode.cpp src\node.cpp src\splitter.cpp src\treedata.cpp src\mtrand.cpp src\datadefs.cpp src\math.cpp src\utils.cpp

cl /EHsc /O2 /analyze /Febin\rf-ace-predict-win32.exe src\rf_ace_predict.cpp src\statistics.cpp src\progress.cpp src\stochasticforest.cpp src\rootnode.cpp src\node.cpp src\splitter.cpp src\treedata.cpp src\mtrand.cpp src\datadefs.cpp src\math.cpp src\utils.cpp

del *.obj

