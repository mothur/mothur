  //
  //  decisiontree.h
  //  rrf-fs-prototype
  //
  //  Created by Abu Zaher Faridee on 5/28/12.
  //  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
  //

#ifndef rrf_fs_prototype_decisiontree_h
#define rrf_fs_prototype_decisiontree_h

#include <iostream>
#include "macros.h"

using namespace std;

class DecisionTree{
public:
  explicit DecisionTree(const vector<TrainingSet>& baseSamples): baseSamples(baseSamples), bootstrappedSamplesSize(baseSamples.size()){
//    DEBUGMSG_VAR(bootstrappedSamplesSize);
    
//    printSamples(baseSamples);
    createBootstrappedSamples();
  }
  
  void createBootstrappedSamples(){
  }
  
    // function that is useful for debugging
  void printSamples(vector<TrainingSet> samples){
    for (unsigned i = 0; i < samples.size(); i++) {
      vector<int> otuCounts = samples[i].getOtuCounts();
      int outputClassId = samples[i].getOutputClassId();
      string outputClass = samples[i].getOutputClass();
      for (unsigned j = 0; j < otuCounts.size(); j++) {
        cout << otuCounts[j] << " ";
      }
      cout << outputClass << " " << outputClassId << endl;
    }
  }
  
private:  
  vector<TrainingSet> baseSamples;
  vector<TrainingSet> bootstrappedTrainingSamples;
  vector<TrainingSet> testSamples;
  
  int bootstrappedSamplesSize;
};


#endif
