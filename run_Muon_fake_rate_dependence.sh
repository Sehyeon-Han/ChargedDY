#!/bin/bash

INPUTFILE="$1"
OUTPUTFILE="$2"

echo "=================================================="
echo "Input file  : ${INPUTFILE}"
echo "Output file : ${OUTPUTFILE}"
echo "Host        : $(hostname)"
echo "Start time  : $(date)"
echo "Working dir : $(pwd)"
echo "=================================================="

source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh

root -l -b -q \
    "Muon_fake_rate_dependence.cc(\"${INPUTFILE}\", \"${OUTPUTFILE}\")"

STATUS=$?

echo "=================================================="
echo "End time    : $(date)"
echo "Exit status : ${STATUS}"
echo "=================================================="

if [ ${STATUS} -ne 0 ]; then
    echo "ERROR: ROOT macro failed."
    exit ${STATUS}
fi

if [ ! -f "${OUTPUTFILE}" ]; then
    echo "ERROR: Output ROOT file was not created."
    exit 1
fi

exit 0
