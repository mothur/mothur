  //
  //  decisiontree.hpp
  //  rrf-fs-prototype
  //
  //  Created by Abu Zaher Faridee on 5/28/12.
  //  Copyright (c) 2012 Schloss Lab. All rights reserved.
  //

#ifndef RF_DECISIONTREE_HPP
#define RF_DECISIONTREE_HPP

#include "macros.h"
#include "rftreenode.hpp"
#include "abstractdecisiontree.hpp"

/***********************************************************************/

struct VariableRankDescendingSorter {
  bool operator() (const pair<int, int>& firstPair, const pair<int, int>& secondPair){
      return firstPair.second > secondPair.second;
  }
};
struct VariableRankDescendingSorterDouble {
    bool operator() (const pair<int, double>& firstPair, const pair<int, double>& secondPair){
        return firstPair.second > secondPair.second;
    }
};
/***********************************************************************/

class DecisionTree: public AbstractDecisionTree{
    
    friend class RandomForest;
    
public:
    
    DecisionTree(vector< vector<int> > baseDataSet,
                 vector<int> globalDiscardedFeatureIndices,
                 OptimumFeatureSubsetSelector optimumFeatureSubsetSelector,
                 string treeSplitCriterion,
                 float featureStandardDeviationThreshold);
    
    virtual ~DecisionTree(){ deleteTreeNodesRecursively(rootNode); }
    
    int calcTreeVariableImportanceAndError(int& numCorrect, double& treeErrorRate);
    int evaluateSample(vector<int> testSample);
    int calcTreeErrorRate(int& numCorrect, double& treeErrorRate);
    
    void randomlyShuffleAttribute(const vector< vector<int> >& samples,
                                  const int featureIndex,
                                  const int prevFeatureIndex,
                                  vector< vector<int> >& shuffledSample);
    
    void purgeDataSetsFromTree() { purgeTreeNodesDataRecursively(rootNode); }
    int purgeTreeNodesDataRecursively(RFTreeNode* treeNode);
    
    void pruneTree(double pruneAggressiveness);
    void pruneRecursively(RFTreeNode* treeNode, double pruneAggressiveness);
    void updateMisclassificationCountRecursively(RFTreeNode* treeNode, vector<int> testSample);
    void updateOutputClassOfNode(RFTreeNode* treeNode);
    
    
private:
    
    void buildDecisionTree();
    int splitRecursively(RFTreeNode* rootNode);
    int findAndUpdateBestFeatureToSplitOn(RFTreeNode* node);
    vector<int> selectFeatureSubsetRandomly(vector<int> globalDiscardedFeatureIndices, vector<int> localDiscardedFeatureIndices);
    int printTree(RFTreeNode* treeNode, string caption);
    void deleteTreeNodesRecursively(RFTreeNode* treeNode);
    
    vector<int> variableImportanceList;
    map<int, int> outOfBagEstimates;
  
    float featureStandardDeviationThreshold;
};

#endif
