  //
  //  dataset.h
  //  rrf-fs-prototype
  //
  //  Created by Abu Zaher Faridee on 5/27/12.
  //  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
  //

#ifndef rrf_fs_prototype_dataset_h
#define rrf_fs_prototype_dataset_h

#include "trainingset.h"

class DataSet{
public:
  explicit DataSet(const vector< vector<string> > sharedFileContent, const vector< vector<string> > designFileContent): 
  sharedFileContent(sharedFileContent),
  designFileContent(designFileContent){
    
    createTrainingSets();
  }
  
//    // copy constructor
//  explicit DataSet(const DataSet& dataSet): 
//    featureLabels(dataSet.featureLabels),
//    numOTUs(dataSet.numOTUs),
//    numTrainingSets(dataSet.numTrainingSets),
//    trainingSets(dataSet.trainingSets),
//    sharedFileContent(dataSet.sharedFileContent),
//    designFileContent(dataSet.designFileContent){
//  }
  
  void createTrainingSets(){
    for(unsigned i = 3; i < sharedFileContent[0].size(); i++){
      featureLabels.push_back(sharedFileContent[0][i]);
    }
    
    string numOUTStr = sharedFileContent[1][2];
    numberOfTotalFeatures = atoi(numOUTStr.c_str());
    
    numberOfTrainingSets = sharedFileContent.size() - 1;
    
    
    for (unsigned i = 1; i < sharedFileContent.size(); i++) {
      
      vector<int> tempOtuCounts;
      for(unsigned j = 3; j < sharedFileContent[i].size(); j++){
        int count = atoi(sharedFileContent[i][j].c_str());
        tempOtuCounts.push_back(count);
      }
      string tempOutputClass(designFileContent[i-1][1].c_str());
      TrainingSet tempTrainingSet(tempOtuCounts, tempOutputClass);
      tempOtuCounts.clear();
        // copy constructor is being called here
      trainingSets.push_back(tempTrainingSet);
    }
    
    createUniqIdForTrainignSets();
//    printTrainingSets();
//    alignTrainingSets();
  }
  
  void printTrainingSets(){
    for(unsigned i = 0; i < trainingSets.size(); i++){
      TrainingSet trainingSet = trainingSets[i];
      for (unsigned j = 0; j < trainingSet.getOtuCounts().size(); j++) {
        cout << trainingSet.getOtuCounts()[j] << " ";
      }
      cout << trainingSet.getOutputClass() << " " << trainingSet.getOutputClassId() << endl;
    }
  }
  
  vector<TrainingSet>& getTrainingSets(){
    return trainingSets;
  }
  
  int getNumberOfTotalFeatures(){
    return numberOfTotalFeatures;
  }

  
private:  
  void createUniqIdForTrainignSets(){
    vector<string> uniqOutputStrings;
    for (unsigned i = 0; i < trainingSets.size(); i++) {
      string outputClass = trainingSets[i].getOutputClass();
      bool found = false;
      for (unsigned j = 0; j < uniqOutputStrings.size(); j++) {
        if (outputClass == uniqOutputStrings[j]){
          found = true;
        }
      }
      if (found == false){
        uniqOutputStrings.push_back(outputClass);
      }
    }
    
    for (unsigned i = 0; i < trainingSets.size(); i++) {
      string outputClass = trainingSets[i].getOutputClass();
      for (unsigned j = 0; j < uniqOutputStrings.size(); j++) {
        if (outputClass == uniqOutputStrings[j]){
          trainingSets[i].setOutputClassId(j);
        }
      }
    }
  }
  
    // some training sets has more data than the rest, the ending cells are 
    // missing in some training sets, so we need to pad those data
  void alignTrainingSets(){
    int maxOtuCounts = -1;
    for (unsigned i = 0; i < trainingSets.size(); i++) {
      int currentOtuCount = trainingSets[i].getOtuCounts().size();
      if (currentOtuCount > maxOtuCounts){ maxOtuCounts = currentOtuCount; }
    }
    
    for (unsigned i = 0; i < trainingSets.size(); i++) {
      int currentOtuCount = trainingSets[i].getOtuCounts().size();
      if (currentOtuCount < maxOtuCounts){
//        cout << "need to add paadding" << endl;
      }
    }
  }
    
  vector<string> featureLabels; 
  unsigned numberOfTotalFeatures; // this is number of OTUs
  unsigned numberOfTrainingSets;
  vector<TrainingSet> trainingSets;
  
  vector< vector<string> > sharedFileContent;
  vector< vector<string> > designFileContent;

};


#endif
