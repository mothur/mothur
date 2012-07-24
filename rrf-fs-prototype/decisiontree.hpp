  //
  //  decisiontree.hpp
  //  rrf-fs-prototype
  //
  //  Created by Abu Zaher Faridee on 5/28/12.
  //  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
  //

#ifndef rrf_fs_prototype_decisiontree_hpp
#define rrf_fs_prototype_decisiontree_hpp

#include <vector>
#include <map>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <limits>
#include "macros.h"
#include "treenode.hpp"
#include "abstractdecisiontree.hpp"

using namespace std;

class DecisionTree: public AbstractDecisionTree{
  
public:
  
  DecisionTree(vector< vector<int> > baseDataSet,
               vector<int> globalDiscardedFeatureIndices,
               OptimumFeatureSubsetSelector optimumFeatureSubsetSelector,
               string treeSplitCriterion)
  
  : AbstractDecisionTree(baseDataSet,
                         globalDiscardedFeatureIndices,
                         optimumFeatureSubsetSelector,
                         treeSplitCriterion),
  
  variableImportanceList(numFeatures, 0){
    
      // TODO self.outOfBagEstimates = {}
    
//		self.buildDecisionTree()
    
    createBootStrappedSamples();
    buildDecisionTree();
  }
  
  ~DecisionTree(){
    if (rootNode != NULL){ delete rootNode; }
  }
  
protected:
  
private:
  
  void buildDecisionTree(){
    int generation = 0;
    rootNode = new TreeNode(bootstrappedTrainingSamples, globalDiscardedFeatureIndices, numFeatures, numSamples, numOutputClasses, generation);
    splitRecursively(rootNode);
#ifdef DEBUG_MODE
//    printTree(rootNode, "root");
#endif
  }
  
  void splitRecursively(TreeNode* rootNode){
  }
  
    // TODO: printTree() needs a check if correct
  void printTree(TreeNode* treeNode, string caption){
    
    string tabs = "";
    for (unsigned i = 0; i < treeNode->getGeneration(); i++) { tabs += "\t"; }
    
    if (treeNode->checkIsLeaf() == false){
      cout << tabs << caption << " [ gen: " << treeNode->getGeneration() << " ] ( " << treeNode->getSplitFeatureValue() << " < X" << treeNode->getSplitFeatureIndex() << " )" << endl;
      printTree(treeNode->getLeftChildNode(), "leftChild");
      printTree(treeNode->getRightChildNode(), "rightChild");
    }else {
      cout << tabs << caption << " [ gen: " << treeNode->getGeneration() + " ] ( classified to: " << treeNode->getOutputClass() << ", samples: " << treeNode->getNumSamples() << " )";
    }
    
  }
  
  vector<int> variableImportanceList;
};

#endif
