#!/usr/bin/bash
echo "start"
X509_USER_PROXY=/afs/cern.ch/user/g/gdamolin/private/x509up_u151129
CMSSW=/afs/cern.ch/user/g/gdamolin/CMSSW_12_4_1_patch1/src
usage() { echo "Usage: $0 [-i <input file> ] [-o <outpath>] [-x xsec] [-l <lumi>] [-d <Data bool>] [-s <Systematcis bool>] [-p <PhysicsProcess string>] [-t test int]" 1>&2; exit 1; }
while getopts "i:o:x:l:d:s:p:t:" opt; do
    case "$opt" in
        i) INFILE=$OPTARG
            ;;
        o) OUTPATH=$OPTARG
            ;;
        p) PRO=$OPTARG
            ;;
        x) XSEC=$OPTARG
            ;;
        l) LUMI=$OPTARG
            ;;
        d) DATA=$OPTARG
            ;;
        s) SYSTS=$OPTARG
            ;;
        t) TEST=$OPTARG
	    ;;
        *)
        usage
        ;;
    esac
done

echo "Process is ${PRO}, Data is $DATA, SYSTEMATICS ARE $SYSTS"
echo "TEST is $TEST"

EXE="/afs/cern.ch/user/g/gdamolin/DYMUMU/DYMuMu/Systs/test.exe"
outdir="/afs/cern.ch/user/g/gdamolin/DYMUMU/SUBDIR"
filename=$INFILE
filestring=$(echo $filename | sed 's|\(^.*/\)\([a-z,A-Z,0-9,-]*\).root$|\2|')
ofilename=${outdir}/$filestring"_MA.root"
echo "ofilename $ofilename"

cd $CMSSW
eval `scram r -sh`
export X509_USER_PROXY=$X509_USER_PROXY
cd -

echo "CMSSW_BASE is now set to $CMSSW_BASE"
echo "PROXY is now set to $X509_USER_PROXY"
echo "executing script as"
echo "${EXE} $filename $ofilename ${XSEC} ${LUMI} ${DATA} ${SYSTS} ${PRO}"
${EXE} $filename $ofilename ${XSEC} ${LUMI} ${DATA} ${SYSTS} ${PRO}|| {
    echo "${EXE} failed with file ${filename}, removing intermediate file" 1>&2; 
    if [[ -f $ofilename ]]; then
        rm $ofilename
    fi
    exit 1
}

echo "$ofilename in ${OUTPATH}/${filestring}_MA.root "
mv $ofilename ${OUTPATH}/${filestring}"_MA.root"


