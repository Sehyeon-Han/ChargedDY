#!/bin/bash

SAMPLE="$1"
INPUTFILE="$2"
MACRO="$3"
OUTPUTFILE="$4"

echo "=================================================="
echo "Sample      : ${SAMPLE}"
echo "Input       : ${INPUTFILE}"
echo "Macro       : ${MACRO}"
echo "Cluster     : ${CLUSTER}"
echo "Process     : ${PROCESS}"
echo "Host        : $(hostname)"
echo "Working dir : $(pwd)"
echo "Start time  : $(date)"
echo "=================================================="


# ============================================================
# ROOT 환경 설정
# 기존 개별 run 파일에서 사용하던 ROOT 설정 경로와 같아야 함
# ============================================================

source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh

if ! command -v root >/dev/null 2>&1; then
    echo "ERROR: ROOT command를 찾을 수 없음"
    exit 10
fi

echo "ROOT version: $(root-config --version)"


# ============================================================
# 입력 검사
# ============================================================

if [ -z "${SAMPLE}" ]; then
    echo "ERROR: SAMPLE이 비어 있음"
    exit 11
fi

if [ -z "${INPUTFILE}" ]; then
    echo "ERROR: INPUTFILE이 비어 있음"
    exit 12
fi

if [ -z "${MACRO}" ]; then
    echo "ERROR: MACRO가 비어 있음"
    exit 13
fi

if [ ! -f "${MACRO}" ]; then
    echo "ERROR: Macro 파일이 없음: ${MACRO}"
    ls -al
    exit 14
fi


# ============================================================
# 출력 파일 이름
# ============================================================

echo "Output file : ${OUTPUTFILE}"


# ============================================================
# ROOT 실행
# ============================================================

echo
echo "[ROOT 실행 시작]"

root -l -b -q "${MACRO}+(\"${INPUTFILE}\",\"${OUTPUTFILE}\")"

ROOT_STATUS=$?

echo "[ROOT 종료 코드: ${ROOT_STATUS}]"

if [ ${ROOT_STATUS} -ne 0 ]; then
    echo "ERROR: ROOT 실행 실패"
    exit ${ROOT_STATUS}
fi


# ============================================================
# 출력 파일 존재 여부 검사
# ============================================================

if [ ! -f "${OUTPUTFILE}" ]; then
    echo "ERROR: 출력 ROOT 파일이 생성되지 않음"
    echo "Expected: ${OUTPUTFILE}"
    ls -al
    exit 20
fi

if [ ! -s "${OUTPUTFILE}" ]; then
    echo "ERROR: 출력 ROOT 파일 크기가 0임"
    rm -f "${OUTPUTFILE}"
    exit 21
fi


# ============================================================
# ROOT 파일 정상 여부 검사
# ============================================================

root -l -b <<EOF

TFile *f = TFile::Open("${OUTPUTFILE}");

if (!f) {
    std::cerr << "ERROR: TFile::Open failed" << std::endl;
    gSystem->Exit(22);
}

if (f->IsZombie()) {
    std::cerr << "ERROR: Zombie ROOT file" << std::endl;
    f->Close();
    delete f;
    gSystem->Exit(23);
}

if (f->TestBit(TFile::kRecovered)) {
    std::cerr << "ERROR: Recovered/corrupted ROOT file" << std::endl;
    f->Close();
    delete f;
    gSystem->Exit(24);
}

std::cout << "ROOT file validation passed" << std::endl;
std::cout << "File size: " << f->GetSize() << " bytes" << std::endl;

f->Close();
delete f;

gSystem->Exit(0);

EOF

CHECK_STATUS=$?

if [ ${CHECK_STATUS} -ne 0 ]; then
    echo "ERROR: ROOT 파일 검사 실패"
    rm -f "${OUTPUTFILE}"
    exit ${CHECK_STATUS}
fi


echo
echo "=================================================="
echo "Job completed successfully"
echo "Output: ${OUTPUTFILE}"
echo "End time: $(date)"
echo "=================================================="

exit 0
