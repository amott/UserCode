#! /bin/csh
# usage: source setup.csh ~/<path>/CMSSW_3_5_6

set CMSSW_release=$1

cd $CMSSW_release
echo setting environment from CMSSW release $CMSSW_release ...
cmsenv 
cd -
#cd TMVA
#echo setting TMVA...
#source setup.csh
#cd ..
echo Done.

