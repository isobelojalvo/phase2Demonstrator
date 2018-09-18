Please note the following can be put into a checkout.sh script and then run using `bash checkout.sh`

```
#!/bin/bash
scramv1 project CMSSW CMSSW_10_1_7
cd CMSSW_10_1_7/src
eval `scramv1 runtime -sh`
git cms-init
git checkout CMSSW_10_1_X
git cms-merge-topic -u isobelojalvo:calocluster-v2
git cms-addpkg L1Trigger/L1TCommon
cd L1Trigger
git clone -b mtd-devel https://github.com/isobelojalvo/phase2Demonstrator.git
cd ../
nohup scramv1 b -j 8 > logBuild &
```