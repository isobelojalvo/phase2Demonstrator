#include "L1Trigger/phase2Demonstrator/interface/triggerGeometryTools.hh"

triggerGeometryTools::triggerGeometryTools(){};



float triggerGeometryTools::getRecoPhi(int iphi){
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
    float temp = towerPhiMap[iphi-1];
    
    return temp;
};

float triggerGeometryTools::getRecoEta(int ieta, short zside){
    //from tower eta -28 to 28
  float towerEtaMap[57]=   { 
      -2.913, //-2.739, 
      -2.565, -2.391, 2.217, 	//switch to smaller trigger towers here    
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

};



int triggerGeometryTools::TPGEtaRange(int ieta)
{
  int iEta = 0;
  // So here, -28 becomes 0.  -1 be comes 27.  +1 becomes 28. +28 becomes 55.
  // And we have mapped [-28, -1], [1, 28] onto [0, 55]   
  if(ieta < 0)
    iEta = ieta + 28;
  else if(ieta > 0)
    iEta = ieta + 27;
  if(ieta==0){
    std::cout<<"Error! ieta is 0, ieta: "<<ieta<<" iEta "<<iEta<<std::endl;
    exit(0);
  }
  return iEta;
}

int triggerGeometryTools::convertGenEta(double inputEta) 
  {
    const double tpgEtaValues[27] = {
      0.087,      
      0.174, // HB and inner HE bins are 0.348 wide
      0.261,
      0.348,
      0.522,
      0.609,
      0.696,
      0.783,
      0.870,
      0.957,
      1.044,
      1.131,
      1.218,
      1.305,
      1.392,
      1.479,
      1.566,
      1.653,
      1.74,
      1.848,
      1.956, // Last two HE bins are 0.432 and 0.828 wide
      2.064,
      2.172,
      2.379,
      2.586,
      2.793,
      3
      //IGNORING HF
      //3.250, // HF bins are 0.5 wide
      //3.750,
      //4.250,
      //4.750
    };

    for (int n=1; n<29; n++){
      //std::cout<<"inputEta "<<inputEta<< " n "<< n <<" tpgEtaValues[n-1] "<< tpgEtaValues[n-1] << " abs(inputEta)<tpgEtaValues[n-1]"<<std::endl;
      if (std::fabs(inputEta)<tpgEtaValues[n-1]) {
	if(inputEta>0){ return n + 20;}
	else{ return n;}
	break;
      }
    }
    std::cout<<"OUT OF BOUNDS!!!!  inputeta: "<<inputEta<<std::endl;
    return -9;
  }

  //-pi < phi <= +pi,
int triggerGeometryTools::convertGenPhi(double inputPhi)
  {
    double posPhi[36];
    for(int n = 0; n < 36; n++)
      posPhi[n] = (0.087) * n + 0.0435;
    double negPhi[36];
    for(int n = 0; n < 36; n++)
      negPhi[n] = -3.14159 + 0.087 * n - 0.0435;

    //1 to 36 is 0 to pi
    if( 3.1416 > inputPhi && inputPhi >= 0){

      for(int n = 1; n < 36; n++){
	//std::cout<<"inputPhi "<<inputPhi<< " posPhi[n-1] "<< posPhi[n-1] << " n "<<n<<std::endl;
	if(inputPhi <= posPhi[n-1]){
	  int tpgPhi = n;
	  return tpgPhi;
	}
      }
    }

    //37 to 72 is -pi to 0
    else if(-3.1416 < inputPhi && inputPhi < 0){
      for(int n = 1; n < 36; n++)
	if(inputPhi < negPhi[n-1]){
	  int tpgPhi = n + 36;
	  return tpgPhi;
	}
    }
    std::cout<<"OUT OF BOUNDS!!!!  inputphi: "<<inputPhi<<std::endl;
    return -9;
  };

uint32_t triggerGeometryTools::getCrystalIEta(float recoEta){
  
  float upperBoundEta[100] = { 0.0174, 0.0348, 0.0522, 0.0696, 0.087, 0.1044, 0.1218, 0.1392, 0.1566, 0.174, 0.1914, 0.2088, 0.2262, 0.2436, 0.261, 0.2784, 0.2958, 0.3132, 0.3306, 0.348, 0.3654, 0.3828, 0.4002, 0.4176, 0.435, 0.4524, 0.4698, 0.4872, 0.5046, 0.522, 0.5394, 0.5568, 0.5742, 0.5916, 0.609, 0.6264, 0.6438, 0.6612, 0.6786, 0.696, 0.7134, 0.7308, 0.7482, 0.7656, 0.783, 0.8004, 0.8178, 0.8352, 0.8526, 0.87, 0.8874, 0.9048, 0.9222, 0.9396, 0.957, 0.9744, 0.9918, 1.0092, 1.0266, 1.044, 1.0614, 1.0788, 1.0962, 1.1136, 1.131, 1.1484, 1.1658, 1.1832, 1.2006, 1.218, 1.2354, 1.2528, 1.2702, 1.2876, 1.305, 1.3224, 1.3398, 1.3572, 1.3746, 1.392, 1.4094, 1.4268, 1.4442, 1.4616, 1.479, 1.4964, 1.5138, 1.5312, 1.5486, 1.566, 1.5834, 1.6008, 1.6182, 1.6356, 1.653, 1.6704, 1.6878, 1.7052, 1.7226, 1.74};

  for(int i = 0; i < 100 ; i++){
    if(abs(recoEta)<upperBoundEta[i])
      return i;
  }
  std::cout<<"ERROR out of bounds eta "<<recoEta<<std::endl;

  return -999;
};

uint32_t triggerGeometryTools::getCrystalIPhi(float recoPhi){
  
  float upperBoundPhi[181] = {0.0174, 0.0348, 0.0522, 0.0696, 0.087, 0.1044, 0.1218, 0.1392, 0.1566, 0.174, 0.1914, 0.2088, 0.2262, 0.2436, 0.261, 0.2784, 0.2958, 0.3132, 0.3306, 0.348, 0.3654, 0.3828, 0.4002, 0.4176, 0.435, 0.4524, 0.4698, 0.4872, 0.5046, 0.522, 0.5394, 0.5568, 0.5742, 0.5916, 0.609, 0.6264, 0.6438, 0.6612, 0.6786, 0.696, 0.7134, 0.7308, 0.7482, 0.7656, 0.783, 0.8004, 0.8178, 0.8352, 0.8526, 0.87, 0.8874, 0.9048, 0.9222, 0.9396, 0.957, 0.9744, 0.9918, 1.0092, 1.0266, 1.044, 1.0614, 1.0788, 1.0962, 1.1136, 1.131, 1.1484, 1.1658, 1.1832, 1.2006, 1.218, 1.2354, 1.2528, 1.2702, 1.2876, 1.305, 1.3224, 1.3398, 1.3572, 1.3746, 1.392, 1.4094, 1.4268, 1.4442, 1.4616, 1.479, 1.4964, 1.5138, 1.5312, 1.5486, 1.566, 1.5834, 1.6008, 1.6182, 1.6356, 1.653, 1.6704, 1.6878, 1.7052, 1.7226, 1.74, 1.7574, 1.7748, 1.7922, 1.8096, 1.827, 1.8444, 1.8618, 1.8792, 1.8966, 1.914, 1.9314, 1.9488, 1.9662, 1.9836, 2.001, 2.0184, 2.0358, 2.0532, 2.0706, 2.088, 2.1054, 2.1228, 2.1402, 2.1576, 2.175, 2.1924, 2.2098, 2.2272, 2.2446, 2.262, 2.2794, 2.2968, 2.3142, 2.3316, 2.349, 2.3664, 2.3838, 2.4012, 2.4186, 2.436, 2.4534, 2.4708, 2.4882, 2.5056, 2.523, 2.5404, 2.5578, 2.5752, 2.5926, 2.61, 2.6274, 2.6448, 2.6622, 2.6796, 2.697, 2.7144, 2.7318, 2.7492, 2.7666, 2.784, 2.8014, 2.8188, 2.8362, 2.8536, 2.871, 2.8884, 2.9058, 2.9232, 2.9406, 2.958, 2.9754, 2.9928, 3.0102, 3.0276, 3.045, 3.0624, 3.0798, 3.0972, 3.1146, 3.132, 3.1494};

  for(int i = 0; i < 181; i++){
    if(abs(recoPhi)<upperBoundPhi[i])
      return i;
  }

  std::cout<<"ERROR out of bounds phi "<<recoPhi<<std::endl;
  
  return -999;

}


uint32_t triggerGeometryTools::getIndex(int towerEta, int towerPhi){
  uint32_t index = 0;

  if(towerEta<0){
    index = abs(towerEta+20)*72 + towerPhi;
  }
  
  if(towerEta>0){
    index = abs(towerEta+20)*72 + towerPhi;
  }
  
  return index;

}
//DEFINE_FWK_MODULE(triggerGeometryTools);
