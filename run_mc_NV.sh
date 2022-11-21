#!/usr/bin/bash
echo "start"
X509_USER_PROXY=/afs/cern.ch/user/g/gdamolin/private/x509up_u151129
CMSSW=/afs/cern.ch/user/g/gdamolin/CMSSW_12_4_1_patch1/src
usage() { echo "Usage: $0 [-i <input file> ] [-o <outpath>] [-x xsec] [-l <lumi>] [-s <Signal bool>]" 1>&2; exit 1; }
while getopts "i:o:x:l:s:" opt; do
    case "$opt" in
        i) INFILE=$OPTARG
            ;;
        o) OUTPATH=$OPTARG
            ;;
        x) XSEC=$OPTARG
            ;;
        l) LUMI=$OPTARG
            ;;
        s) SIGNAL=$OPTARG
            ;;
        *)
        usage
        ;;
    esac
done
EXE="/afs/cern.ch/user/g/gdamolin/DYMUMU/DYMuMu/Noveto_MC.exe"
outdir="/afs/cern.ch/user/g/gdamolin/DYMUMU/noveto"
filename=$INFILE
filestring=$(echo $filename | sed 's|\(^.*/\)\([a-z,A-Z,0-9,-]*\).root$|\2|')
ofilename=${outdir}/$filestring"_MA".root

cd $CMSSW
eval `scram r -sh`
export X509_USER_PROXY=$X509_USER_PROXY
cd -

echo "CMSSW_BASE is now set to $CMSSW_BASE"
echo "PROXY is now set to $X509_USER_PROXY"
echo "executing script as"
echo "${EXE} $filename $ofilename ${XSEC} ${LUMI} ${SIGNAL}"
${EXE} $filename $ofilename ${XSEC} ${LUMI} ${SIGNAL} || {
    echo "${EXE} failed with file ${filename}, removing intermediate file" 1>&2; 
    if [[ -f $ofilename ]]; then
        rm $ofilename
    fi
    exit 1
}
mv $ofilename ${OUTPATH}/${filestring}_MA.root
