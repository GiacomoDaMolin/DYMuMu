# DYMuMu
Study of DY->MuMu at CMS
Hitchhiker’s guide to the galaxy
0) create afs submit directory (where .dag and rescue are placed) and eos (root file storage)
1) create dag files with submit_dag_ROOT.py
python3 submit_dag_ROOT.py -i /afs/cern.ch/user/g/gdamolin/DYMUMU/SUBDIR/ -o /eos/user/g/gdamolin/DYMUMU/ -j /afs/cern.ch/user/g/gdamolin/DYMUMU/DYMuMu/samples_MUDT.json -e /afs/cern.ch/user/g/gdamolin/DYMUMU/DYMuMu/run_data.sh -p /afs/cern.ch/user/g/gdamolin/private/x509up_u151129
or
python3 submit_dag_ROOT.py -m -i /afs/cern.ch/user/g/gdamolin/DYMUMU/SUBDIR/ -o /eos/user/g/gdamolin/DYMUMU/ -j /afs/cern.ch/user/g/gdamolin/DYMUMU/DYMuMu/samples_mcv2.json -e /afs/cern.ch/user/g/gdamolin/DYMUMU/DYMuMu/run_mc.sh -p /afs/cern.ch/user/g/gdamolin/private/x509up_u151129

2) use
for dir in $(find ./ -maxdepth 1 -mindepth 1 -type d); do condor_submit_dag $dir/*.dag; done
to launch dag first time

3) if some process fails, use rescue.py
python3 rescue.py -s /afs/cern.ch/user/g/gdamolin/DYMUMU/SUBDIR/ -o /eos/user/g/gdamolin/DYMUMU/ -j /afs/cern.ch/user/g/gdamolin/DYMUMU/DYMuMu/samples_mc_v2.json (-r)
