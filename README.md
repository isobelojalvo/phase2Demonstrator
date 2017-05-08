# phase2Demonstrator
L1TriggerPhase2Demonstrator

```
scramv1 project CMSSW CMSSW_9_1_0_pre2
cd CMSSW_9_1_0_pre2/src
eval `scramv1 runtime -sh`
cmsenv
git cms-merge-topic isobelojalvo:CMSSW_9_0_X_L1PF_devel
git cms-merge-topic skinnari:Tracklet_91X
cd L1Trigger
git clone git@github.com:isobelojalvo/phase2Demonstrator.git
cd ../
nohup scramv1 b -j 8
```