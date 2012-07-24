#!/bin/sh

prefix=$1

echo "SUBMITTING DATA..."
python cmst3_submit_manyfilesperjob_eeAnalysis.py Data7TeV dataset_eg_Sep3rdReReco 10 VecbosApp 4 8nh $prefix 0
echo "DONE WITH DATA."

