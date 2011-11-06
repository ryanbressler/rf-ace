mkdir bin
cl /EHsc /Ox /Febin\rf_ace_win64.exe src\rf_ace.cpp src\stochasticforest.cpp src\rootnode.cpp src\node.cpp src\treedata.cpp src\mtrand.cpp src\datadefs.cpp src\partitionsequence.cpp
del *.obj