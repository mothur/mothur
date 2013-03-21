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
  bool operator() (vector<int> first, vector<int> second){ return first[1] > second[1]; }
};
struct VariableRankDescendingSorterDouble {
    bool operator() (vector<double> first, vector<double> second){ return first[1] > second[1]; }
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
    vector< vector<int> > randomlyShuffleAttribute(vector< vector<int> > samples, int featureIndex);  
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
