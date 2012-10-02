SOURCEFILES: rf_ace_R.cpp progress.cpp statistics.cpp math.cpp stochasticforest.cpp rootnode.cpp node.cpp treedata.cpp datadefs.cpp utils.cpp distributions.cpp

all: rfacer

rfacer: $(SOURCEFILES)
	$(CPP) -shared $(PKG_LIBS) $(PKG_CPPFLAGS) $(SOURCEFILES) -pthread -o rfacer.so
