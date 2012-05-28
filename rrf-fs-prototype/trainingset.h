  //
  //  trainingset.h
  //  rrf-fs-prototype
  //
  //  Created by Abu Zaher Faridee on 5/27/12.
  //  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
  //

#ifndef rrf_fs_prototype_trainingset_h
#define rrf_fs_prototype_trainingset_h

class TrainingSet{
public:
    // explicit constructor
  explicit TrainingSet(vector<int> otuCounts, string outputClass): 
      otuCounts(otuCounts),
      outputClass(outputClass),
      outputClassId(0){
  }
  
    // copy constructor
  TrainingSet(const TrainingSet& trainingSet): 
      otuCounts(trainingSet.otuCounts),
      outputClass(trainingSet.outputClass),
      outputClassId(trainingSet.outputClassId){
  }
  
  vector<int>& getOtuCounts(){
    return otuCounts;
  }
  
  string getOutputClass(){
    return outputClass;
  }
  
  int getOutputClassId(){
    return outputClassId;
  }
  
  void setOutputClass(string outputClass){
    this->outputClass = outputClass;
  }
  
  void setOutputClassId(int outputClassId){
    this->outputClassId = outputClassId;
  }
  
private:
  vector<int> otuCounts;
  string outputClass;
  int outputClassId;
};


#endif
