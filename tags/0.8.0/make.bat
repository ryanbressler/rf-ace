mkdir bin
cl /Ox /Febin\rf_ace.exe src\rf_ace.cpp src\stochasticforest.cpp src\rootnode.cpp src\node.cpp src\treedata.cpp src\mtrand.cpp src\datadefs.cpp
del *.obj