  //
  //  decisiontree.h
  //  rrf-fs-prototype
  //
  //  Created by Abu Zaher Faridee on 5/28/12.
  //  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
  //

#ifndef rrf_fs_prototype_decisiontree_h
#define rrf_fs_prototype_decisiontree_h

class DecisionTree{
public:
  explicit DecisionTree(DataSet dataset):
  dataSet(dataSet),
  bootstrappedSamplesSize(100){
    dataSet.printTrainingSets();
      //      cout << bootstrappedSamplesSize << endl;
  }
  
  void createBootstrappedSamples(){
  }
private:
  DataSet dataSet;
  vector<TrainingSet> bootstrappedTrainingSamples;
  vector<TrainingSet> testSamples;
  
  int bootstrappedSamplesSize;
};


#endif
