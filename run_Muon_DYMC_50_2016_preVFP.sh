#!/bin/bash

source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh

FILE=$1 
IDX=$2
OUTFILE=Muon_DYMC_50_2016_preVFP_${IDX}.root 

echo "Running on $FILE" 

root -l "Muon_DYMC_50_2016_preVFP.cc(\"$FILE\",\"$OUTFILE\")" 