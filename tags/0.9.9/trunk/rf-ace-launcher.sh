#!/bin/sh
# Generic launcher for RF-ACE
#  
#  This script is used primarily with GenePattern, allowing us to specify
#  all options explicitly, with defaults provided by GenePattern's web 
#  interface. In effect, this makes all parameters positional.
#
#  The GenePattern command line for this is:
#   sh <libdir>rf-ace-launcher.sh <libdir> <input> <target> <output> \
#    <RF_ntrees> <RF_mtry> <RF_nodesize> <RF_nperms> <RF_pthreshold> \
#    <GBT_ntrees> <GBT_maxleaves> <GBT_shrinkage> <GBT_samplesize> \
#    "<RF_enable><RF_optimize><GBT_enable><GBT_optimize><fmask>"

export PATH=$1:$PATH

chmod a+x $1/rf_ace
echo "Running: \
rf_ace --input=$2 --target=$3 --output=$4 --RF_ntrees=$5 \
--RF_mtry=$6 --RF_nodesize=$7 --RF_nperms=$8 --RF_pthreshold=$9 \
--GBT_ntrees=${10} --GBT_maxleaves=${11} --GBT_shrinkage=${12} \
--GBT_samplesize=${13} ${14}"

rf_ace --input=$2 --target=$3 --output=$4 --RF_ntrees=$5 \
--RF_mtry=$6 --RF_nodesize=$7 --RF_nperms=$8 --RF_pthreshold=$9 \
--GBT_ntrees=${10} --GBT_maxleaves=${11} --GBT_shrinkage=${12} \
--GBT_samplesize=${13} ${14}
