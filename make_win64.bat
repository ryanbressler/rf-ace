mkdir bin
cl /EHsc /Ox /analyze /Febin\rf-ace-filter-win64.exe src\rf_ace_filter.cpp src\statistics.cpp src\progress.cpp src\stochasticforest.cpp src\rootnode.cpp src\node.cpp src\splitter.cpp src\treedata.cpp src\mtrand.cpp src\datadefs.cpp src\utils.cpp src\math.cpp src\partitionsequence.cpp
del *.obj