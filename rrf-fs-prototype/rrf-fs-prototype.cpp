  //
  //  main.cpp
  //  rrf-fs-prototype
  //
  //  Created by Abu Zaher Faridee on 5/21/12.
  //  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
  //
  //  Note: This is a very basic implementation of Regularized Random Forest 
  //  Algorithm for Feature Selection. After the prototype has been done, 
  //  we'd need to integrate this into mothur as a seperate module.
  //

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "sharedanddesignfilereader.h"
#include "dataset.h"

class TrainingSet;
class DataSet;

class RegularizedRandomForest{
public:
  RegularizedRandomForest(string sharedFilePath, string designFilePath) : 
    sharedFileReader(sharedFilePath), 
    designFileReader(designFilePath),
    dataSet(sharedFileReader.getFileContent(), designFileReader.getFileContent()){
          
//    sharedFileReader.printFileContent();    
//    designFileReader.printFileContent();
    
    dataSet.createTrainingSets();
//    dataSet.printTrainingSets();

  }
  
  void createBootstrappedSamples(){
  }
private:
  
  SharedAndDesignFileReader sharedFileReader;
  SharedAndDesignFileReader designFileReader;
  
  DataSet dataSet;
};


using namespace std;

int main(int argc, const char * argv[]){
  
  string sharedFilePath = "final.an.0.03.subsample.0.03.pick.shared";
  string designFilePath = "mouse.sex_time.design";
  
  RegularizedRandomForest regularizedRandomForest(sharedFilePath, designFilePath);
  
  return 0;
}

