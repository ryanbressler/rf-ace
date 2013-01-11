mkdir bin

SetEnv.cmd /x86 /Release

cl /EHsc /O2 /analyze /DNOTHREADS /Febin\rf-ace-win32.exe src\murmurhash3.cpp src\rf_ace.cpp src\statistics.cpp src\distributions.cpp src\progress.cpp src\stochasticforest.cpp src\rootnode.cpp src\node.cpp src\treedata.cpp src\datadefs.cpp src\math.cpp src\utils.cpp src\reader.cpp src\feature.cpp

del *.obj

