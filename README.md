#!/bin/bash
scramv1 project CMSSW CMSSW_10_1_7
cd CMSSW_10_1_7/src
eval `scramv1 runtime -sh`
git cms-init
git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git
git fetch cms-l1t-offline phase2-l1t-integration-CMSSW_10_1_7
git cms-merge-topic -u cms-l1t-offline:l1t-phase2-v2.16.24
git cms-merge-topic -u isobelojalvo:l1caloclusters
cd L1Trigger
git clone -b mtd-devel https://github.com/isobelojalvo/phase2Demonstrator.git
cd ../
nohup scramv1 b -j 8