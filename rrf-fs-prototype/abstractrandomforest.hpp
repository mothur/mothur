//
//  abstractrandomforest.hpp
//  rrf-fs-prototype
//
//  Created by Abu Zaher Faridee on 7/20/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef rrf_fs_prototype_abstractrandomforest_hpp
#define rrf_fs_prototype_abstractrandomforest_hpp

#include <vector>
#include <map>
#include <numeric>
#include <algorithm>
#include <cmath>
#include "macros.h"
#include "abstractdecisiontree.hpp"

#define DEBUG_MODE

using namespace std;


class AbstractRandomForest{
public:
    // intialization with vectors
  AbstractRandomForest(const std::vector < std::vector<int> > dataSet, 
                       const int numDecisionTrees, 
                       const std::string treeSplitCriterion = "informationGain")
  : dataSet(dataSet), 
  numDecisionTrees(numDecisionTrees),
  numSamples(dataSet.size()),
  numFeatures(dataSet[0].size() - 1),
  globalDiscardedFeatureIndices(getGlobalDiscardedFeatureIndices()),
  globalVariableImportanceList(numFeatures),
  treeSplitCriterion(treeSplitCriterion){
      // TODO: double check if the implemenatation of 'globalOutOfBagEstimates' is correct
  }
  
    // intialization with 2d const array
  AbstractRandomForest(int** dataSetAs2dArray, const int rows, const int columns, const int numDecisionTrees, const std::string treeSplitCriterion="informationGain"){
  }
    
  ~AbstractRandomForest(){
  }
  
protected:
  
  vector<int> getGlobalDiscardedFeatureIndices(){
      // calculate feature vectors
    vector< vector<int> > featureVectors(numFeatures);
    for (int i = 0; i < numSamples; i++) {
      for (int j = 0; j < numFeatures; j++) {
        featureVectors[j].push_back(dataSet[i][j]);
      }
    }
    
    vector<int> globalDiscardedFeatureIndices;
    
    for (int i = 0; i < featureVectors.size(); i++) {
      double standardDeviation = getStandardDeviation(featureVectors[i]);
      if (standardDeviation <= 0){ globalDiscardedFeatureIndices.push_back(i); }
    }
    
#ifdef DEBUG_MODE
    cout << "number of global discarded features: "<< globalDiscardedFeatureIndices.size() << endl;
    cout << "total features: " << featureVectors.size() << endl;
#endif
    return globalDiscardedFeatureIndices;
  }
  
private:
  int numDecisionTrees;
  int numSamples;
  int numFeatures;
  vector< vector<int> > dataSet;
  vector<int> globalDiscardedFeatureIndices;
  vector<double> globalVariableImportanceList;
  string treeSplitCriterion;
  map<int, vector<int> > globalOutOfBagEstimates;
  vector<AbstractDecisionTree> decisionTrees;
};
#endif
