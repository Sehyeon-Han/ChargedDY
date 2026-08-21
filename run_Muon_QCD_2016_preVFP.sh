#!/bin/bash

INPUTFILE=$1
OUTPUTFILE=$2

echo "=================================================="
echo "Input       : ${INPUTFILE}"
echo "Output      : ${OUTPUTFILE}"
echo "Host        : $(hostname)"
echo "Working dir : $(pwd)"
echo "Start time  : $(date)"
echo "=================================================="

source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh

echo "ROOT version:"
root-config --version

ls -lh

root -l -b -q "Muon_QCD_2016_preVFP.cc(\"${INPUTFILE}\",\"${OUTPUTFILE}\",\"real_mu_rate.root\",\"fake_mu_normalized_results.root\")"

ROOT_STATUS=$?

if [ ${ROOT_STATUS} -ne 0 ]; then
    echo "ERROR: ROOT command failed with status ${ROOT_STATUS}"
    exit ${ROOT_STATUS}
fi

if [ ! -f "${OUTPUTFILE}" ]; then
    echo "ERROR: Output file was not created"
    exit 1
fi

echo "=================================================="
echo "Output file:"
ls -lh "${OUTPUTFILE}"
echo "End time: $(date)"
echo "=================================================="

exit 0