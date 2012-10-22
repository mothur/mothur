//
//  macros.h
//  rrf-fs-prototype
//
//  Created by Abu Zaher Faridee on 5/28/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#ifndef rrf_fs_prototype_macros_h
#define rrf_fs_prototype_macros_h

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
