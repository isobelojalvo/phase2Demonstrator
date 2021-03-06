#ifndef triggerGeometryTools_hh
#define triggerGeometryTools_hh

//#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
//#include "DataFormats/L1TrackTrigger/interface/TTPixelTrack.h"
//#include "DataFormats/L1TrackTrigger/interface/L1TkPrimaryVertex.h"

#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"

#include <vector>
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <TLorentzVector.h>

class triggerGeometryTools{
public:

  triggerGeometryTools();
  //input is ieta
  float getRecoEta(int ieta, short zside);
  
  //input is iphi
  float getRecoPhi(int iphi);  
  
  int TPGEtaRange(int ieta);

  int convertGenEta(double inputEta);

  int convertGenPhi(double inputPhi);
  
  float SumTLorentzPt(std::vector<TLorentzVector> inputVector){
    float sum = 0;
    for(auto element : inputVector){
      sum += element.Pt();
    }
    return(sum);
  }

  uint32_t getCrystalIEta(float recoEta);  
  uint32_t getCrystalIPhi(float recoPhi);

  uint32_t getIndex(int towerEta, int towerPhi);

  ~triggerGeometryTools(){}

protected:
  

};


#endif
