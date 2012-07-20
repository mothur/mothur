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
#include "regularizedrandomforest.h"
#include "Datasets/inpatient.final.an.0.03.subsample.avg.matrix.h"

class TrainingSet;
class DataSet;


using namespace std;

int main(int argc, const char * argv[]){
  
    // call srand only once in the program
  srand(time(NULL));
  
//  const string sharedFilePath = "final.an.0.03.subsample.0.03.pick.shared";
//  const string designFilePath = "mouse.sex_time.design";
//  const int numberOfDecisionTrees = 1000;
  
//  RegularizedRandomForest regularizedRandomForest(sharedFilePath, designFilePath, numberOfDecisionTrees);
  
    // test of Dataset
//  for (int i = 0; i < rows; i++) {
//    for (int j = 0; j < columns; j++) {
//      cout << dataSet[i][j] << " ";
//    }
//    cout << endl;
//  }
      
  return 0;
}

