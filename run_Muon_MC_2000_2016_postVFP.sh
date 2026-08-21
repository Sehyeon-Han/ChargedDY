#!/bin/bash

source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh

FILE=$1 
IDX=$2
OUTFILE=Muon_MC_2000_2016_postVFP_${IDX}.root 

echo "Running on $FILE" 

root -l "Muon_MC_2000_2016_postVFP.cc(\"$FILE\",\"$OUTFILE\")" 