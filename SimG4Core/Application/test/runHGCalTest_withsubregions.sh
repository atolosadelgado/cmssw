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
cmsRun CMSSW_15_1_0/src/SimG4Core/Application/test/HGCalTest_withsubregions.py
