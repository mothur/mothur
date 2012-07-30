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
#include <limits>

#include "Datasets/inpatient.final.an.0.03.subsample.avg.matrix.h"
#include "randomforest.hpp"
#include "decisiontree.hpp"
#include "treenode.hpp"


using namespace std;

//bool comparator(vector<int> first, vector<int> second){ return first[1] < second[1]; }

int main(int argc, const char * argv[]){
  
// call srand only once in the program
  srand((unsigned int)time(NULL));
  
//  const string sharedFilePath = "final.an.0.03.subsample.0.03.pick.shared";
//  const string designFilePath = "mouse.sex_time.design";
//  const int numberOfDecisionTrees = 1000;
      
  /* create 2d vector from the array */
  vector< vector<int> > dataSet(numRows, vector<int>(numColumns, 0));
  for (int i = 0; i < numRows; i++) { for (int j = 0; j < numColumns; j++) { dataSet[i][j] = inpatientDataSet[i][j]; } }
  
  int numDecisionTrees = 1;
  RandomForest randomForest(dataSet, numDecisionTrees, "informationGain");
      
  randomForest.populateDecisionTrees();
  randomForest.calcForrestErrorRate();
  randomForest.calcForrestVariableImportance();

  
//   depending on the first element of the row, the 2d matrix is sorted
//   the whole row is shuffled not just the first element of the row
//  vector< vector<int> > a(3, vector<int>(3, 0));
//  a[0][0] = 3; a[0][1] = 18; a[0][2] = 7;
//  a[1][0] = 2; a[1][1] = 12; a[1][2] = 1;
//  a[2][0] = 14; a[2][1] = 1; a[2][2] = -1;
//
//  PRINT_VAR(a);  
//  sort(a.begin(), a.end());
//  PRINT_VAR(a);
//  sort(a.begin(), a.end(), comparator);
//  PRINT_VAR(a);
//  sort(a.begin(), a.end(), obj);
//  PRINT_VAR(a);
  
}