#!/usr/bin/sh

package=$1

tar -czf $package src/*.*pp bin/rf-ace* test/*.*pp test/matlab/*.m test/python/*.py test_*by*.tsv test_*by*.arff test_*by*.afm test_predictor.sf testdata.tsv Makefile make_win32.bat make_win64.bat doxy.cfg rf_ace_batch.py rf-ace-launcher.sh 