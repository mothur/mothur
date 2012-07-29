//
//  randomforest.hpp
//  rrf-fs-prototype
//
//  Created by Abu Zaher Faridee on 7/20/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef rrf_fs_prototype_randomforest_hpp
#define rrf_fs_prototype_randomforest_hpp

#include <vector>
#include <map>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <limits>
#include "macros.h"
#include "abstractrandomforest.hpp"
#include "decisiontree.hpp"

using namespace std;

class RandomForest: public AbstractRandomForest {
  
public:
  
    // DONE
  RandomForest(const vector <vector<int> > dataSet,
               const int numDecisionTrees,
               const string treeSplitCriterion = "informationGain")
  
  : AbstractRandomForest(dataSet, numDecisionTrees, treeSplitCriterion) {
    
  }
  
  virtual ~RandomForest() {
    for (vector<AbstractDecisionTree*>::iterator it = decisionTrees.begin(); it != decisionTrees.end(); it++) {
        // we know that this is decision tree, so we can do a dynamic_case<DecisionTree*> here
      DecisionTree* decisionTree = dynamic_cast<DecisionTree*>(*it);
        // calling the destructor by deleting
      delete decisionTree;
    }
  }
  
    // DONE
  void calcForrestErrorRate() {
    
#ifdef DEBUG_LEVEL_2
    DEBUGMSG_FUNC;
#endif
    
    int numCorrect = 0;
    for (map<int, vector<int> >::iterator it = globalOutOfBagEstimates.begin(); it != globalOutOfBagEstimates.end(); it++) {
      int indexOfSample = it->first;
      vector<int> predictedOutComes = it->second;
      vector<int>::iterator maxPredictedOutComeIterator = max_element(predictedOutComes.begin(), predictedOutComes.end());
      int majorityVotedOutcome = (int)(maxPredictedOutComeIterator - predictedOutComes.begin());
      int realOutcome = dataSet[indexOfSample][numFeatures];
      
      if (majorityVotedOutcome == realOutcome) { numCorrect++; }
    }
    
      // TODO: save or return forrestErrorRate for future use;
    double forrestErrorRate = 1 - ((double)numCorrect / (double)globalOutOfBagEstimates.size());
    
#ifdef DEBUG_LEVEL_2
    PRINT_VAR(globalOutOfBagEstimates.size());
    PRINT_VAR(numCorrect);
    PRINT_VAR(forrestErrorRate);
#endif
    
  }
  
    // DONE
  void calcForrestVariableImportance() {
    
#ifdef DEBUG_MODE
    DEBUGMSG_FUNC;
#endif
    
      // TODO: need to add try/catch operators to fix this
      // follow the link: http://en.wikipedia.org/wiki/Dynamic_cast
    for (unsigned i = 0; i < decisionTrees.size(); i++) {
      DecisionTree* decisionTree = dynamic_cast<DecisionTree*>(decisionTrees[i]);
      
      for (unsigned j = 0; j < numFeatures; j++) {
        globalVariableImportanceList[j] += (double)decisionTree->variableImportanceList[j];
      }
    }
    
    for (unsigned i = 0;  i < numFeatures; i++) {
      globalVariableImportanceList[i] /= (double)numDecisionTrees;
    }
    
    vector< vector<int> > globalVariableRanks;
    for (unsigned i = 0; i < globalVariableImportanceList.size(); i++) {
      if (globalVariableImportanceList[i] > 0) {
        vector<int> globalVariableRank(2, 0);
        globalVariableRank[0] = i; globalVariableRank[1] = globalVariableImportanceList[i];
        globalVariableRanks.push_back(globalVariableRank);
      }
    }
    
    VariableRankDescendingSorter variableRankDescendingSorter;
    sort(globalVariableRanks.begin(), globalVariableRanks.end(), variableRankDescendingSorter);
    
#ifdef DEBUG_MODE
    PRINT_VAR(globalVariableRanks);
#endif
    
  }
  
    // DONE
  void populateDecisionTrees() {
    
#ifdef DEBUG_MODE
    DEBUGMSG_FUNC;
#endif
    
    
    for (unsigned i = 0; i < numDecisionTrees; i++) {
      cout << "Creating " << i << " (th) Decision tree" << endl;
        // TODO: need to first fix if we are going to use pointer based system or anything else
      DecisionTree* decisionTree = new DecisionTree(dataSet, globalDiscardedFeatureIndices, OptimumFeatureSubsetSelector("log2"), treeSplitCriterion);
      decisionTree->calcTreeVariableImportanceAndError();
      updateGlobalOutOfBagEstimates(decisionTree);
      decisionTree->purgeDataSetsFromTree();
      decisionTrees.push_back(decisionTree);
    }
    
#ifdef DEBUG_MODE
    PRINT_VAR(globalOutOfBagEstimates);
#endif
    
  }
  
    // TODO: need to finalize bettween reference and pointer for DecisionTree [partially solved]
    // TODO: make this pure virtual in superclass
    // DONE
  void updateGlobalOutOfBagEstimates(DecisionTree* decisionTree) {
    for (map<int, int>::iterator it = decisionTree->outOfBagEstimates.begin(); it != decisionTree->outOfBagEstimates.end(); it++) {
      int indexOfSample = it->first;
      int predictedOutcomeOfSample = it->second;
      
      if (globalOutOfBagEstimates.count(indexOfSample) == 0) {
        globalOutOfBagEstimates[indexOfSample] = vector<int>(decisionTree->numOutputClasses, 0);
      };
      
      globalOutOfBagEstimates[indexOfSample][predictedOutcomeOfSample] += 1;
    }
  }

  
protected:
  
private:
  
};

#endif
