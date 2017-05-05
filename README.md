# phase2Demonstrator
L1TriggerPhase2Demonstrator

```
scramv1 project CMSSW CMSSW_9_1_0_pre2
cd CMSSW_9_1_0_pre2/src
eval `scramv1 runtime -sh`
cmsenv
git cms-merge-topic skinnari:Tracklet_91X

nohup scramv1 b -j 8
```