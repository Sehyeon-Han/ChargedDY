#!/bin/bash

OUTPUT="all_inputs_preVFP.txt"

rm -f "${OUTPUT}"
touch "${OUTPUT}"

add_sample() {
    SAMPLE="$1"
    INPUT_LIST="$2"
    MACRO="$3"

    if [ ! -f "${INPUT_LIST}" ]; then
        echo "[WARNING] input list 없음: ${INPUT_LIST}"
        return
    fi

    if [ ! -f "${MACRO}" ]; then
        echo "[WARNING] C++ macro 없음: ${MACRO}"
        return
    fi

    COUNT=0

    while IFS= read -r INPUTFILE || [ -n "${INPUTFILE}" ]; do

        # 앞뒤 공백 제거
        INPUTFILE=$(echo "${INPUTFILE}" | xargs)

        # 빈 줄 제외
        [ -z "${INPUTFILE}" ] && continue

        # 주석 줄 제외
        [[ "${INPUTFILE}" == \#* ]] && continue

        echo "${SAMPLE} ${INPUTFILE} ${MACRO}" >> "${OUTPUT}"

        COUNT=$((COUNT + 1))

    done < "${INPUT_LIST}"

    echo "[ADD] ${SAMPLE}: ${COUNT} jobs"
}


# ============================================================
# W -> mu nu signal
# ============================================================

add_sample \
    "MC_200" \
    "input_WToMuNu_M-200_2016_preVFP.txt" \
    "Muon_MC_200_2016_preVFP.cc"

add_sample \
    "MC_500" \
    "input_WToMuNu_M-500_2016_preVFP.txt" \
    "Muon_MC_500_2016_preVFP.cc"

add_sample \
    "MC_1000" \
    "input_WToMuNu_M-1000_2016_preVFP.txt" \
    "Muon_MC_1000_2016_preVFP.cc"

add_sample \
    "MC_2000" \
    "input_WToMuNu_M-2000_2016_preVFP.txt" \
    "Muon_MC_2000_2016_preVFP.cc"


# ============================================================
# W -> tau nu
# ============================================================

add_sample \
    "WToTauNu_M-200" \
    "input_WToTauNu_M-200_2016_preVFP.txt" \
    "Muon_WToTauNu_M_200_2016_preVFP.cc"

add_sample \
    "WToTauNu_M-500" \
    "input_WToTauNu_M-500_2016_preVFP.txt" \
    "Muon_WToTauNu_M_500_2016_preVFP.cc"

add_sample \
    "WToTauNu_M-1000" \
    "input_WToTauNu_M-1000_2016_preVFP.txt" \
    "Muon_WToTauNu_M_1000_2016_preVFP.cc"

add_sample \
    "WToTauNu_M-2000" \
    "input_WToTauNu_M-2000_2016_preVFP.txt" \
    "Muon_WToTauNu_M_2000_2016_preVFP.cc"


# ============================================================
# TTbar
# ============================================================

add_sample \
    "TTTo2L2Nu" \
    "input_TTTo2L2Nu_2016_preVFP.txt" \
    "Muon_TTTo2L2Nu_2016_preVFP.cc"

add_sample \
    "TTToSemiLeptonic" \
    "input_TTToSemiLeptonic_2016_preVFP.txt" \
    "Muon_TTToSemiLeptonic_2016_preVFP.cc"


# ============================================================
# Single Top
# ============================================================

add_sample \
    "ST_s_channel" \
    "input_ST_s-channel_4f_leptonDecays_2016_preVFP.txt" \
    "Muon_ST_s_channel_4f_leptonDecays_2016_preVFP.cc"

add_sample \
    "ST_t_antitop" \
    "input_ST_t-channel_anitop_4f_InclusiveDecays_2016_preVFP.txt" \
    "Muon_ST_t_channel_antitop_4f_InclusiveDecays_2016_preVFP.cc"

add_sample \
    "ST_t_top" \
    "input_ST_t-channel_top_4f_InclusiveDecays_2016_preVFP.txt" \
    "Muon_ST_t_channel_top_4f_InclusiveDecays_2016_preVFP.cc"

add_sample \
    "ST_tW_antitop" \
    "input_ST_tW_antitop_5f_inclusiveDecays_2016_preVFP.txt" \
    "Muon_ST_tW_antitop_5f_inclusiveDecays_2016_preVFP.cc"

add_sample \
    "ST_tW_top" \
    "input_ST_tW_top_5f_inclusiveDecays_2016_preVFP.txt" \
    "Muon_ST_tW_top_5f_inclusiveDecays_2016_preVFP.cc"


# ============================================================
# WW
# ============================================================

add_sample \
    "WWTo1L1Nu2Q" \
    "input_WWTo1L1Nu2Q_2016_preVFP.txt" \
    "Muon_WWTo1L1Nu2Q_2016_preVFP.cc"

add_sample \
    "WWTo2L2Nu" \
    "input_WWTo2L2Nu_2016_preVFP.txt" \
    "Muon_WWTo2L2Nu_2016_preVFP.cc"

add_sample \
    "WWTo4Q" \
    "input_WWTo4Q_4f_2016_preVFP.txt" \
    "Muon_WWTo4Q_4f_2016_preVFP.cc"


# ============================================================
# WZ
# ============================================================

add_sample \
    "WZTo1L1Nu2Q" \
    "input_WZTo1L1Nu2Q_4f_2016_preVFP.txt" \
    "Muon_WZTo1L1Nu2Q_4f_2016_preVFP.cc"

add_sample \
    "WZTo1L3Nu" \
    "input_WZTo1L3Nu_4f_2016_preVFP.txt" \
    "Muon_WZTo1L3Nu_4f_2016_preVFP.cc"

add_sample \
    "WZTo2Q2Nu" \
    "input_WZTo2Q2Nu_4f_2016_preVFP.txt" \
    "Muon_WZTo2Q2Nu_4f_2016_preVFP.cc"

add_sample \
    "WZTo3LNu" \
    "input_WZTo3LNu_2016_preVFP.txt" \
    "Muon_WZTo3LNu_2016_preVFP.cc"

# ============================================================
# ZZ
# ============================================================

add_sample \
    "ZZ" \
    "input_ZZ_2016_preVFP.txt" \
    "Muon_ZZ_2016_preVFP.cc"


echo
echo "============================================"
echo "생성 파일: ${OUTPUT}"
echo "총 job 개수: $(wc -l < "${OUTPUT}")"
echo "============================================"