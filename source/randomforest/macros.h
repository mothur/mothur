//
//  macros.h
//  rrf-fs-prototype
//
//  Created by Abu Zaher Faridee on 5/28/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#ifndef RF_MACROS_H
#define RF_MACROS_H

#include "mothurout.h" 

/***********************************************************************/
class OptimumFeatureSubsetSelector{
public:
  OptimumFeatureSubsetSelector(string selectionType = "log2"): selectionType(selectionType){
  }
  
  int getOptimumFeatureSubsetSize(int numFeatures){

    if (selectionType == "log2"){ return (int)ceil(log2(numFeatures)); }
    else if (selectionType == "squareRoot"){ return (int)ceil(sqrt(numFeatures)); } 
    return -1;
  }
private:
  string selectionType;
};

/***********************************************************************/
  
#endif
