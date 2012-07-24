#!/bin/csh

# Get the TMVA version 4.0.3 from SVN
svn co https://tmva.svn.sourceforge.net/svnroot/tmva/tags/V04-00-03/TMVA TMVA

# compile TMVA
cd TMVA
source setup.csh
cd src
make 

# Come back to the original directory
cd ../../
echo "*******************************************************************"
echo "TMVA version 4.0.3 installed in :$PWD/TMVA"
echo "*******************************************************************"

