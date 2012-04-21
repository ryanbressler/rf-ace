mkdir bin

SetEnv.cmd /x64 /Release

cl /EHsc /O2 /analyze /Febin\rf-ace-win64.exe src\rf_ace.cpp src\statistics.cpp src\progress.cpp src\stochasticforest.cpp src\rootnode.cpp src\node.cpp src\treedata.cpp src\mtrand.cpp src\datadefs.cpp src\math.cpp src\gamma.cpp src\utils.cpp

del *.obj

