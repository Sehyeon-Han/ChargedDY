#!/bin/bash

source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh

FILE=$1 
IDX=$2
OUTFILE=Muon_WWTo2L2Nu_2016_preVFP_${IDX}.root 

echo "Running on $FILE" 

root -l "Muon_WWTo2L2Nu_2016_preVFP.cc(\"$FILE\",\"$OUTFILE\")" 