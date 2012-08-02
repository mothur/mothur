  //
  //  trainingset.h
  //  rrf-fs-prototype
  //
  //  Created by Abu Zaher Faridee on 5/27/12.
  //  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
  //

#ifndef rrf_fs_prototype_trainingset_h
#define rrf_fs_prototype_trainingset_h

#include "macros.h"

class TrainingSet{
public:
    // explicit constructor
  explicit TrainingSet(const vector<int> otuCounts, const string outputClass): 
      otuCounts(otuCounts),
      outputClass(outputClass),
      outputClassId(0){
  }
  
    // copy constructor
  TrainingSet(const TrainingSet& trainingSet): 
    otuCounts(trainingSet.otuCounts),
    outputClass(trainingSet.outputClass),
    outputClassId(trainingSet.outputClassId){
        
//    cout << "TrainingSet copy constructor is being called" << endl;
//    for (int i = 0; i < otuCounts.size(); i++) {
//      cout << otuCounts[i] << " ";
//    }
//    cout << outputClass << " " << outputClassId << endl;
        
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
  
  void setOutputClass(const string outputClass){
    this->outputClass = outputClass;
  }
  
  void setOutputClassId(const int outputClassId){
    this->outputClassId = outputClassId;
  }

  friend ostream& operator <<(ostream& os, TrainingSet& trainingSet);
  
private:
  vector<int> otuCounts;
  string outputClass;
  int outputClassId;
};

ostream& operator <<(ostream& os, TrainingSet& trainingSet){
  os << "{ ";
  NAME_VALUE_PAIR(trainingSet.otuCounts, os);
  os << ",\n";
  NAME_VALUE_PAIR(trainingSet.outputClass, os);
  os << ",\n";
  NAME_VALUE_PAIR(trainingSet.outputClassId, os);
  os << " }";
  return os;
}

ostream& operator <<(ostream& os, vector<TrainingSet>& trainingSets){
  os << "[ ";
  for (int i = 0; i < trainingSets.size(); i++) {
    os << trainingSets[i] << "\n\n";
  }
  os << "]";
  return os;
}


#endif
