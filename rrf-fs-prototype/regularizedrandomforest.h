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
  explicit RegularizedRandomForest(const string sharedFilePath, const string designFilePath, const int numberOfDecisionTrees) : 
    sharedFileReader(sharedFilePath), 
    designFileReader(designFilePath),
    dataSet(sharedFileReader.getFileContent(), designFileReader.getFileContent()),
    numberOfDecisionTrees(numberOfDecisionTrees){
      
//    sharedFileReader.printFileContent();    
//    designFileReader.printFileContent();    
//    dataSet.printTrainingSets();
      
    vector<TrainingSet> baseSamples = dataSet.getTrainingSets();
      
//    cout << baseSamples << endl;
      
    for (unsigned i = 0; i < numberOfDecisionTrees; i++) {
      DecisionTree decisionTree(baseSamples);
      decisionTrees.push_back(decisionTree);
    }
  
  }
  
private:
  
  SharedAndDesignFileReader sharedFileReader;
  SharedAndDesignFileReader designFileReader;
  
  DataSet dataSet;
  vector<DecisionTree> decisionTrees;
  
  int numberOfDecisionTrees;
};


#endif
