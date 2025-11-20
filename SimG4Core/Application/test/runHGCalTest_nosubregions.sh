#!/bin/bash
set -e

#set condor ids
CLUSTERID=$1
PROCID=$2

if [ -z "$CLUSTERID" ] || [ -z "$PROCID" ]; then
    echo "ERROR: This script requires 2 arguments: <ClusterId> <ProcId>"
    exit 1
fi

git config --global user.github atolosadelgado
localdir=$PWD
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SITECONFIG_PATH=/cvmfs/cms.cern.ch/SITECONF/T0_CH_CERN
cmsrel CMSSW_15_1_0
cd CMSSW_15_1_0/src
cmsenv

git cms-addpkg SimG4Core/Application
git cms-addpkg SimG4Core/PhysicsLists
git fetch my-cmssw
git checkout my-cmssw/eprofile-allevents
scram b -j4
cd $localdir


cmsRun CMSSW_15_1_0/src/SimG4Core/Application/test/HGCalTest_nosubregions.py

BASE_OUTPUT_DIR="/afs/cern.ch/user/a/atolosad/work/HGCal/condor/test_nosubregions"
JOB_OUTPUT_DIR="${BASE_OUTPUT_DIR}/${CLUSTERID}.${PROCID}"
mv *.root $JOB_OUTPUT_DIR
