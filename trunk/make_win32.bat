mkdir bin
cl /EHsc /Ox /Febin\rf_ace_win32.exe src\rf_ace.cpp src\statistics.cpp src\progress.cpp src\stochasticforest.cpp src\rootnode.cpp src\node.cpp src\splitter.cpp src\treedata.cpp src\mtrand.cpp src\datadefs.cpp src\utils.cpp src\partitionsequence.cpp
del *.obj