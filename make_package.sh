#!/usr/bin/sh

package=$1

tar -czf $package src/*.*pp bin/rf-ace* test/ test_*.* testdata.tsv Makefile make_win32.bat make_win64.bat doxy.cfg rf_ace_batch.py rf-ace-launcher.sh 
