//
//  abstractdecisiontree.hpp
//  rrf-fs-prototype
//
//  Created by Abu Zaher Faridee on 7/22/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef rrf_fs_prototype_abstractdecisiontree_hpp
#define rrf_fs_prototype_abstractdecisiontree_hpp

#include <vector>
#include <map>
#include <numeric>
#include <algorithm>
#include <cmath>
#include "macros.h"

#define DEBUG_MODE

using namespace std;

class AbstractDecisionTree{
public:
  AbstractDecisionTree(vector<vector<int> >baseDataSet, 
                       vector<int> globalDiscardedFeatureIndices, 
                       OptimumFeatureSubsetSelector optimumFeatureSubsetSelector, 
                       string treeSplitCriterion)
  : baseDataSet(baseDataSet),
  numSamples(baseDataSet.size()),
  numFeatures(baseDataSet[0].size()),
  numOutputClasses(0),
    // TODO: self.rootNode = None
  globalDiscardedFeatureIndices(globalDiscardedFeatureIndices),
  optimumFeatureSubsetSize(optimumFeatureSubsetSelector.getOptimumFeatureSubsetSize(numFeatures)),
  treeSplitCriterion(treeSplitCriterion){
    
    for (int i = 0;  i < numSamples; i++) {
      int outcome = baseDataSet[i][numFeatures - 1];
      vector<int>::iterator it = find(classes.begin(), classes.end(), outcome);
      if (it == classes.end()){       // find() will return classes.end() if the element is not found
        classes.push_back(outcome);
        numOutputClasses++;
      }
    }    
  }
  
  ~AbstractDecisionTree(){
  }
protected:
private:
  vector< vector<int> > baseDataSet;
  int numSamples;
  int numFeatures;
  int numOutputClasses;
  vector<int> classes;
  vector< vector<int> > bootstrappedTrainingSamples;
  vector<int> bootstrappedTrainingSampleIndices;
  vector< vector<int> > bootstrappedTestSamples;
  vector<int> bootstrappedTestSampleIndices;
    // self.rootNode = None
  vector<int> globalDiscardedFeatureIndices;
  int optimumFeatureSubsetSize; // fix this
  string treeSplitCriterion;
};

#endif
