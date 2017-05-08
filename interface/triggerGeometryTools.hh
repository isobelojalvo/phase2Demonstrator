#ifndef triggerGeometryTools_hh
#define triggerGeometryTools_hh

#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"

#include <vector>
#include <iostream>

#include "L1Trigger/L1TCaloLayer1/src/UCTRegion.hh"
#include "L1Trigger/L1TCaloLayer1/src/UCTGeometry.hh"


/* */
float upperBoundCrystalEta[86] = {  0.0174, 0.0348, 0.0522, 0.0696, 0.087, 0.1044, 0.1218, 0.1392, 0.1566, 0.174, 0.1914, 0.2088, 0.2262, 0.2436, 0.261, 0.2784, 0.2958, 0.3132, 0.3306, 0.348, 0.3654, 0.3828, 0.4002, 0.4176, 0.435, 0.4524, 0.4698, 0.4872, 0.5046, 0.522, 0.5394, 0.5568, 0.5742, 0.5916, 0.609, 0.6264, 0.6438, 0.6612, 0.6786, 0.696, 0.7134, 0.7308, 0.7482, 0.7656, 0.783, 0.8004, 0.8178, 0.8352, 0.8526, 0.87, 0.8874, 0.9048, 0.9222, 0.9396, 0.957, 0.9744, 0.9918, 1.0092, 1.0266, 1.044, 1.0614, 1.0788, 1.0962, 1.1136, 1.131, 1.1484, 1.1658, 1.1832, 1.2006, 1.218, 1.2354, 1.2528, 1.2702, 1.2876, 1.305, 1.3224, 1.3398, 1.3572, 1.3746, 1.392, 1.4094, 1.4268, 1.4442, 1.4616, 1.479, 1.4964
};


//from tower eta -28 to 28
float towerEtaMap[57]=   { 
  -2.913, //-2.739, 
  -2.565, -2.391, 2.217, //switch to smaller trigger towers here    
  -2.0445, -1.9575, -1.8705, -1.7835, -1.6965, 
  -1.6095, -1.5225, -1.4355, -1.3485, -1.2615, 
  -1.1745, -1.0875, -1.0005, -0.9135, -0.8265, 
  -0.7395, -0.6525, -0.5655, -0.4785, -0.3915, 
  -0.3045, -0.2175, -0.1305, -0.0435, 0.0435, 
  0.1305, 0.2175, 0.3045, 0.3915, 0.4785, 
  0.5655, 0.6525, 0.7395, 0.8265, 0.9135, 
  1.0005, 1.0875, 1.1745, 1.2615, 1.3485, 
  1.4355, 1.5225, 1.6095, 1.6965, 1.7835, 
  1.8705, 1.9575, 2.0445, 2.217, 2.391, 
  2.565, //2.739,
    2.913 
};

float towerPhiMap[72]=                        
  {0.044, 0.131, 0.218, 0.305, 0.393, 0.480, 0.567, 0.654, 0.742, 0.829, 0.916, 1.004, 1.091, 1.178, 1.265, 1.353, 1.440, 1.527, 1.614, 1.702, 1.789, 1.876, 1.963, 2.051, 2.138, 2.225, 2.313, 2.400, 2.487, 2.574, 2.662, 2.749, 2.836, 2.923, 3.011, 3.098,
   -3.098, -3.011, -2.923, -2.836, -2.749, -2.662, -2.574, -2.487, -2.400, -2.313, -2.225, -2.138, -2.051, -1.963, -1.876, -1.789, -1.702, -1.614, -1.527, -1.440, -1.353, -1.265, -1.178, -1.091, -1.004, -0.916, -0.829, -0.742, -0.654, -0.567, -0.480, -0.393, -0.305, -0.218, -0.131, -0.044};
/*{-0.1305, -0.0435,  0.0435, 0.1308, 0.2178, 0.3048, 
   0.3918,  0.4788, 0.5658, 0.6528, 0.7398, 0.8268, 
   0.9138, 1.0008, 1.0878, 1.1748, 1.2618, 1.3488, 
   1.4358, 1.5228, 1.6098, 1.6968, 1.7838, 1.8708, 
   1.9578, 2.0448, 2.1318, 2.2188, 2.3058, 2.3928, 
   2.4798, 2.5668, 2.6538, 2.7408, 2.8278, 2.9148, 
   3.0018, 3.0888, -3.0885, -3.0015, -2.9145, -2.8275, 
   -2.7405, -2.6535, -2.5665, -2.4795, -2.3925, -2.3055, 
   -2.2185, -2.1315, -2.0445, -1.9575, -1.8705, -1.7835, 
   -1.6965, -1.6095, -1.5225, -1.4355, -1.3485, -1.2615, 
   -1.1745, -1.0875, -1.0005, -0.9135, -0.8265, -0.7395, 
   -0.6525, -0.5655, -0.4785, -0.3915, -0.3045, -0.2175 };*/



//input is iphi
float getRecoPhi(int iphi){
  return towerPhiMap[iphi-1];
}

//input is ieta
float getRecoEta(int ieta, short zside){
  float eta = -999;
  if(ieta<0 || ieta>(28*2)){
    std::cout<<"Error!!! towereta out of bounds in triggerGeometryTools.h "<<std::endl;
    std::cout<<"ieta "<<ieta<<std::endl;
    exit(0);
  }
  
  if(zside == 1)
    eta = towerEtaMap[ieta];
  else if(zside == -1)
    eta = towerEtaMap[ieta];
  else{
    std::cout<<"Error!!! zside out of bounds in triggerGeometryTools.h "<<std::endl;
    std::cout<<"zside "<<zside<<std::endl;
    exit(0);
  }
  return eta;
}



#endif

