#!/bin/bash

# =========================
# 사용자 설정
# =========================

INPUT_LIST="input_WToMuNu_M-200_2016_preVFP.txt"
OUTPUT_DIR="."
OUTPUT_PREFIX="Muon_MC_200_2016_preVFP_"
OUTPUT_SUFFIX=".root"

RETRY_SUBMIT="submit_retry.sub"
MISSING_LIST="missing_jobs.txt"

# 1이면 파일 존재 여부만 검사
# 0이면 ROOT 파일이 정상적으로 열리는지도 검사
CHECK_EXISTENCE_ONLY=0


# =========================
# 초기 설정
# =========================

if [ ! -f "$INPUT_LIST" ]; then
    echo "[ERROR] input list가 없습니다: $INPUT_LIST"
    exit 1
fi

if [ ! -f "$RETRY_SUBMIT" ]; then
    echo "[ERROR] retry submit 파일이 없습니다: $RETRY_SUBMIT"
    exit 1
fi

mkdir -p logs_retry

> "$MISSING_LIST"

TOTAL=$(grep -cv '^[[:space:]]*$' "$INPUT_LIST")

echo "========================================"
echo "전체 예상 job 수: $TOTAL"
echo "검사 시작"
echo "========================================"


# =========================
# output 검사
# =========================

jobidx=0

while IFS= read -r filename; do

    # 빈 줄 무시
    [ -z "$filename" ] && continue

    output_file="${OUTPUT_DIR}/${OUTPUT_PREFIX}${jobidx}${OUTPUT_SUFFIX}"

    retry=0
    reason=""

    # 파일이 없거나 크기가 0인 경우
    if [ ! -s "$output_file" ]; then
        retry=1
        reason="MISSING 또는 EMPTY"

    # ROOT 파일 정상 여부 검사
    elif [ "$CHECK_EXISTENCE_ONLY" -eq 0 ]; then

        root -l -b > /dev/null 2>&1 <<EOF
TFile *f = TFile::Open("$output_file");

if (!f || f->IsZombie() || f->TestBit(TFile::kRecovered)) {
    if (f) f->Close();
    gSystem->Exit(1);
}

f->Close();
gSystem->Exit(0);
EOF

        if [ $? -ne 0 ]; then
            retry=1
            reason="BROKEN ROOT"
        fi
    fi

    if [ "$retry" -eq 1 ]; then
        echo "[RETRY] job $jobidx : $reason"
        echo "$jobidx $filename" >> "$MISSING_LIST"
    else
        echo "[OK] job $jobidx"
    fi

    jobidx=$((jobidx + 1))

done < "$INPUT_LIST"


# =========================
# 결과 및 재제출
# =========================

MISSING_COUNT=$(wc -l < "$MISSING_LIST")

echo
echo "========================================"
echo "누락 또는 손상 job 수: $MISSING_COUNT"
echo "========================================"

if [ "$MISSING_COUNT" -eq 0 ]; then
    echo "모든 output이 정상적으로 생성되었습니다."
    rm -f "$MISSING_LIST"
    exit 0
fi

echo
echo "재제출 대상:"
cat "$MISSING_LIST"

echo
echo "다음 submit 파일로 재제출합니다:"
echo "$RETRY_SUBMIT"
echo

condor_submit "$RETRY_SUBMIT"

if [ $? -eq 0 ]; then
    echo
    echo "누락된 job 재제출 완료"
else
    echo
    echo "[ERROR] condor_submit 실패"
    exit 1
fi