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
  numSamples((int)dataSet.size()),
  numFeatures((int)(dataSet[0].size() - 1)),
  globalDiscardedFeatureIndices(getGlobalDiscardedFeatureIndices()),
  globalVariableImportanceList(numFeatures, 0),
  treeSplitCriterion(treeSplitCriterion){
      // TODO: double check if the implemenatation of 'globalOutOfBagEstimates' is correct
  }
  
    // intialization with 2d const array
//  AbstractRandomForest(int** dataSetAs2dArray, const int rows, const int columns, const int numDecisionTrees, const std::string treeSplitCriterion="informationGain"){
//  }
    
  virtual ~AbstractRandomForest(){
  }
  
  virtual void populateDecisionTrees() = 0;
  virtual void calcForrestErrorRate() = 0;
  virtual void calcForrestVariableImportance() = 0;
  
protected:
  
    // TODO: create a better way of discarding feature
    // currently we just set FEATURE_DISCARD_SD_THRESHOLD to 0 to solved this
    // it can be tuned for better selection
    // also, there might be other factors like Mean or other stuffs
    // same would apply for createLocalDiscardedFeatureList in the TreeNode class
  
    // TODO: Another idea is getting an aggregated discarded feature indices after the run, from combining
    // the local discarded feature indices
    // this would penalize a feature, even if in global space the feature looks quite good
    // the penalization would be averaged, so this woould unlikely to create a local optmina
  vector<int> getGlobalDiscardedFeatureIndices() {
      // calculate feature vectors
    vector< vector<int> > featureVectors(numFeatures, vector<int>(numSamples, 0));
    for (int i = 0; i < numSamples; i++) {
      for (int j = 0; j < numFeatures; j++) { featureVectors[j][i] = dataSet[i][j]; }
    }
    
    vector<int> globalDiscardedFeatureIndices;
    
    for (int i = 0; i < featureVectors.size(); i++) {
      double standardDeviation = getStandardDeviation(featureVectors[i]);
      if (standardDeviation <= FEATURE_DISCARD_SD_THRESHOLD){ globalDiscardedFeatureIndices.push_back(i); }
    }
    
#ifdef DEBUG_MODE
    PRINT_MSG("number of global discarded features: ", globalDiscardedFeatureIndices.size());
    PRINT_MSG("total features:", featureVectors.size());
#endif
    
    return globalDiscardedFeatureIndices;
  }

  int numDecisionTrees;
  int numSamples;
  int numFeatures;
  vector< vector<int> > dataSet;
  vector<int> globalDiscardedFeatureIndices;
  vector<double> globalVariableImportanceList;
  string treeSplitCriterion;
    // This is a map of each feature to outcome count of each classes
    // e.g. 1 => [2 7] means feature 1 has 2 outcome of 0 and 7 outcome of 1
  map<int, vector<int> > globalOutOfBagEstimates;
  
    // TODO: fix this, do we use pointers?
  vector<AbstractDecisionTree*> decisionTrees;
  
private:

};
#endif
