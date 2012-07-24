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
  
    // TODO: splitRecursively() finish implementation
  void splitRecursively(TreeNode* rootNode){
    
#ifdef DEBUG_MODE
    DEBUGMSG_LOCATION;
#endif
    
    if (rootNode->getNumSamples() < 2){      
#ifdef DEBUG_MODE
      DEBUGMSG("Already classified: Case 1");
#endif      
      rootNode->setIsLeaf(true);
      rootNode->setOutputClass(rootNode->getBootstrappedTrainingSamples()[0][rootNode->getNumFeatures()]);
      return;
    }
    
    int classifiedOutputClass;
    bool isAlreadyClassified = checkIfAlreadyClassified(rootNode, classifiedOutputClass);    
    if (isAlreadyClassified == true){
#ifdef DEBUG_MODE
      DEBUGMSG("Already classified: Case 2");
#endif
      rootNode->setIsLeaf(true);
      rootNode->setOutputClass(classifiedOutputClass);
      return;
    }
    
    vector<int> featureSubsetIndices = selectFeatureSubsetRandomly(globalDiscardedFeatureIndices, rootNode->getLocalDiscardedFeatureIndices());
    rootNode->setFeatureSubsetIndices(featureSubsetIndices);
    
#ifdef DEBUG_MODE
    DEBUGMSG_VAR(globalDiscardedFeatureIndices);
    DEBUGMSG_VAR(featureSubsetIndices);
#endif
    
    findAndUpdateBestFeatureToSplitOn(rootNode);
    
    vector< vector<int> > leftChildSamples;
    vector< vector<int> > rightChildSamples;
    getSplitPopulation(rootNode, leftChildSamples, rightChildSamples);
    
      // TODO: need to write code to clear this memory
    TreeNode* leftChildNode = new TreeNode(leftChildSamples, globalDiscardedFeatureIndices, numFeatures, leftChildSamples.size(), numOutputClasses, rootNode->getGeneration() + 1);
    TreeNode* rightChildNode = new TreeNode(rightChildSamples, globalDiscardedFeatureIndices, numFeatures, rightChildSamples.size(), numOutputClasses, rootNode->getGeneration() + 1);
    
    rootNode->setLeftChildNode(leftChildNode);
    leftChildNode->setParentNode(rootNode);
    
    rootNode->setRightChildNode(rightChildNode);
    rightChildNode->setParentNode(rootNode);
    
      // TODO: This recursive split can be parrallelized later
    splitRecursively(leftChildNode);
    splitRecursively(rightChildNode);
  }
  
    // TODO: finish implementation of findAndUpdateBestFeatureToSplitOn()
  void findAndUpdateBestFeatureToSplitOn(TreeNode* rootNode){
  }
  
  vector<int> selectFeatureSubsetRandomly(vector<int> globalDiscardedFeatureIndices, vector<int> localDiscardedFeatureIndices){
#ifdef DEBUG_MODE
    DEBUGMSG_LOCATION;
#endif
    vector<int> featureSubsetIndices;
    
    vector<int> combinedDiscardedFeatureIndices;
    combinedDiscardedFeatureIndices.insert(combinedDiscardedFeatureIndices.end(), globalDiscardedFeatureIndices.begin(), globalDiscardedFeatureIndices.end());
    combinedDiscardedFeatureIndices.insert(combinedDiscardedFeatureIndices.end(), localDiscardedFeatureIndices.begin(), localDiscardedFeatureIndices.end());
    
    sort(combinedDiscardedFeatureIndices.begin(), combinedDiscardedFeatureIndices.end());
    
    int numberOfRemainingSuitableFeatures = numFeatures - combinedDiscardedFeatureIndices.size();
    int currentFeatureSubsetSize = numberOfRemainingSuitableFeatures < optimumFeatureSubsetSize ? numberOfRemainingSuitableFeatures : optimumFeatureSubsetSize;
    
    while (featureSubsetIndices.size() < currentFeatureSubsetSize) {
      int randomIndex = rand() % numFeatures;
      vector<int>::iterator it = find(featureSubsetIndices.begin(), featureSubsetIndices.end(), randomIndex);
      if (it == featureSubsetIndices.end()){    // NOT FOUND
        vector<int>::iterator it2 = find(combinedDiscardedFeatureIndices.begin(), combinedDiscardedFeatureIndices.end(), randomIndex);
        if (it2 == combinedDiscardedFeatureIndices.end()){  // NOT FOUND AGAIN
          featureSubsetIndices.push_back(randomIndex);
        }
      }
    }
    sort(featureSubsetIndices.begin(), featureSubsetIndices.end());
    
    return featureSubsetIndices;
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
