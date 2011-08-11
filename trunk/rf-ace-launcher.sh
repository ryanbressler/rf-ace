#!/bin/sh
# Generic launcher for RF-ACE
#  
#  This script is used primarily with GenePattern, allowing us to specify
#  all options explicitly, with defaults provided by GenePattern's web 
#  interface. In effect, this makes all parameters positional.
#
#  The GenePattern command line for this is:
#   sh <libdir>rf-ace-launcher.sh <libdir> <input> <target> <output> \
#    <ntrees> <mtry> <nodesize> <nperms> <pthreshold> <maxleaves> \
#    <shrinkage> <subsamplesize> <filter_rf> <optimized_rf> <enable_gbt>
#

export PATH=$1:$PATH
flags=
if [ ${13} -eq "1" ]; then flags="--filterRF "; fi
if [ ${14} -eq "1" ]; then flags="$flags--optimizedRF "; fi
if [ ${15} -eq "1" ]; then flags="$flags--enableGBT "; fi

chmod a+x $1/rf_ace
echo "Running: rf_ace --input=$2 --target=$3 --output=$4 --ntrees=$5 --mtry=$6 \
--nodesize=$7 --nperms=$8 --pthreshold=$9 --maxleaves=${10} --shrinkage=${11} \
--subsamplesize=${12} $flags"

rf_ace --input=$2 --target=$3 --output=$4 --ntrees=$5 --mtry=$6 \
--nodesize=$7 --nperms=$8 --pthreshold=$9 --maxleaves=${10} --shrinkage=${11} \
--subsamplesize=${12} $flags
