  //
  //  regularizedrandomforest.h
  //  rrf-fs-prototype
  //
  //  Created by Abu Zaher Faridee on 5/28/12.
  //  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
  //

#ifndef rrf_fs_prototype_regularizedrandomforest_h
#define rrf_fs_prototype_regularizedrandomforest_h

#include "decisiontree.h"

class RegularizedRandomForest{
public:
  explicit RegularizedRandomForest(string sharedFilePath, string designFilePath, int numberOfDecisionTrees) : 
  sharedFileReader(sharedFilePath), 
  designFileReader(designFilePath),
  dataSet(sharedFileReader.getFileContent(), designFileReader.getFileContent()),
  numberOfDecisionTrees(numberOfDecisionTrees){
    
    sharedFileReader.printFileContent();    
    designFileReader.printFileContent();    
    dataSet.printTrainingSets();
    
      //    DecisionTree decisionTree(dataSet);
    
  }
  
private:
  
  SharedAndDesignFileReader sharedFileReader;
  SharedAndDesignFileReader designFileReader;
  
  DataSet dataSet;
  vector<DecisionTree> decisionTress;
  
  int numberOfDecisionTrees;
};


#endif
