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
#include "treenode.hpp"

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
  rootNode(NULL),
  globalDiscardedFeatureIndices(globalDiscardedFeatureIndices),
  optimumFeatureSubsetSize(optimumFeatureSubsetSelector.getOptimumFeatureSubsetSize(numFeatures)),
  treeSplitCriterion(treeSplitCriterion){
    
    for (int i = 0;  i < numSamples; i++) {
      int outcome = baseDataSet[i][numFeatures - 1];
      vector<int>::iterator it = find(outputClasses.begin(), outputClasses.end(), outcome);
      if (it == outputClasses.end()){       // find() will return classes.end() if the element is not found
        outputClasses.push_back(outcome);
        numOutputClasses++;
      }
    }    
  }
  
  ~AbstractDecisionTree(){
    if (rootNode != NULL){ delete rootNode; }
  }
protected:
private:
  vector< vector<int> > baseDataSet;
  int numSamples;
  int numFeatures;
  int numOutputClasses;
  vector<int> outputClasses;
  vector< vector<int> > bootstrappedTrainingSamples;
  vector<int> bootstrappedTrainingSampleIndices;
  vector< vector<int> > bootstrappedTestSamples;
  vector<int> bootstrappedTestSampleIndices;
  
  TreeNode* rootNode;
  vector<int> globalDiscardedFeatureIndices;
  int optimumFeatureSubsetSize;
  string treeSplitCriterion;
};

#endif
