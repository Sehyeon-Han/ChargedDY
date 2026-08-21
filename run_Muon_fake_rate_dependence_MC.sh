#!/bin/bash

set -e

SAMPLE=$1
INPUTFILE=$2
OUTPUTFILE=$3

echo "=================================================="
echo "Sample      : ${SAMPLE}"
echo "Input       : ${INPUTFILE}"
echo "Output      : ${OUTPUTFILE}"
echo "Host        : $(hostname)"
echo "Working dir : $(pwd)"
echo "Start time  : $(date)"
echo "=================================================="

source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh

echo "ROOT version:"
root-config --version

if [ ! -f "Muon_fake_rate_dependence_MC.cc" ]; then
    echo "ERROR: Muon_fake_rate_dependence_MC.cc does not exist."
    exit 1
fi

root -l -b -q \
"Muon_fake_rate_dependence_MC.cc(\"${INPUTFILE}\",\"${OUTPUTFILE}\")"

if [ ! -f "${OUTPUTFILE}" ]; then
    echo "ERROR: Output ROOT file was not created."
    exit 2
fi

echo "=================================================="
echo "Finished successfully"
echo "Output size : $(du -h "${OUTPUTFILE}" | cut -f1)"
echo "End time    : $(date)"
echo "=================================================="
