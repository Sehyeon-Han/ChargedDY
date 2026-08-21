#!/bin/bash

FILE=$1
OUT=$2

source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh

root -l -b -q "sumW.cc(\"$FILE\", \"$OUT\")"