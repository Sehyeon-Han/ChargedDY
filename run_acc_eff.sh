#!/bin/bash
# Usage: run_acc_eff.sh MASS INPUT_FILE CLUSTER_ID PROCESS_ID
set -e
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 MASS INPUT_FILE CLUSTER_ID PROCESS_ID" >&2
    exit 2
fi
mass=$1
input_file=$2
cluster_id=$3
process_id=$4
case "$mass" in
    100|200|500|1000|2000) ;;
    *) echo "Unsupported mass: $mass" >&2; exit 2 ;;
esac
if [[ ! "$cluster_id" =~ ^[0-9]+$ || ! "$process_id" =~ ^[0-9]+$ ]]; then
    echo "Cluster and process IDs must be integers" >&2
    exit 2
fi
source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh
set -u
macro="acc_eff_M${mass}.cc"
out_file="acc_eff_M${mass}_${cluster_id}_${process_id}.root"
for required in "$macro" eff_mu_ID.root eff_mu_ISO.root eff_mu_STRIG.root PUreweight.root; do
    if [ ! -s "$required" ]; then
        echo "Missing or empty required file: $required" >&2
        exit 3
    fi
done
# Escape paths for C++ string literals passed to ROOT.
cpp_input=${input_file//\\/\\\\}
cpp_input=${cpp_input//\"/\\\"}
echo "Input: $input_file"
echo "Output: $out_file"
root -l -b -q "${macro}+(\"${cpp_input}\",\"${out_file}\")"
if [ ! -s "$out_file" ]; then
    echo "Missing or empty output: $out_file" >&2
    exit 4
fi
root -l -b <<EOF_ROOT
TFile *f = TFile::Open("${out_file}");
if (!f || f->IsZombie() || f->TestBit(TFile::kRecovered)) { gSystem->Exit(5); }
if (!f->Get("h_acc_den") || !f->Get("h_acc_num") || !f->Get("h_eff_num")) { gSystem->Exit(6); }
f->Close();
gSystem->Exit(0);
EOF_ROOT
echo "Completed: $out_file"
