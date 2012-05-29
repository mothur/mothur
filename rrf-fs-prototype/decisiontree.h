  //
  //  decisiontree.h
  //  rrf-fs-prototype
  //
  //  Created by Abu Zaher Faridee on 5/28/12.
  //  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
  //

#ifndef rrf_fs_prototype_decisiontree_h
#define rrf_fs_prototype_decisiontree_h

#include <iostream>
#include <cmath>
#include "macros.h"
#include "treenode.hpp"

using namespace std;

class DecisionTree{
public:
  explicit DecisionTree(const vector<TrainingSet>& baseSamples, int numberOfTotalFeatures): 
    baseSamples(baseSamples), 
    bootstrappedSamplesSize(baseSamples.size()),
    numberOfTotalFeatures(numberOfTotalFeatures){
//    DEBUGMSG_VAR(bootstrappedSamplesSize);
    
    createBootstrappedSamples();
    buildDecisionTree();
  }
  
    // need to manually provide a copy constructor so that when addign a DecisionTree to a vector, we get desired results
  DecisionTree(const DecisionTree& decisionTree) : 
    baseSamples(decisionTree.baseSamples),
    bootstrappedTrainingSamples(decisionTree.bootstrappedTrainingSamples),
    bootStrammpedTestSamples(decisionTree.bootStrammpedTestSamples),
    bootstrappedSamplesSize(decisionTree.bootstrappedSamplesSize),
    numberOfTotalFeatures(decisionTree.numberOfTotalFeatures){
  }
  
  void buildDecisionTree(){
    TreeNode* rootTreeNode = new TreeNode(bootstrappedTrainingSamples);
    splitRecursively(rootTreeNode);
    delete rootTreeNode;
  }
  
  void splitRecursively(TreeNode* rootTreeNode){
    if(rootTreeNode->isLeaf == false){
      vector<int> featureSubsetIndices = selectFeatureSubsetRandomly();
    }
  }

    // TODO: we need to modify this later for randomized random forest algo's implementations
  vector<int> selectFeatureSubsetRandomly(){
    vector<int> featureSubsetIndices;
    int optimumFeatureSubsetSize = getOptimumFeatureSubsetSize();
//    DEBUGMSG_VAR(optimumFeatureSubsetSize);
    
    vector<bool> isSelected(numberOfTotalFeatures, false);
//    DEBUGMSG_VAR(isSelected);
    
    while (true) {
        // TODO: be careful with this random index calculation function, there is a chance it will not
        // output the hightest number say (527) ever.
      double randomIndex = (int)((double)rand() / (double)RAND_MAX * numberOfTotalFeatures);
      isSelected[randomIndex] = true;
      
      int count = 0;
      for (unsigned i = 0; i < numberOfTotalFeatures; i++) {
        if (isSelected[i] == true){ 
          count++; 
        }
      }
      if (count == optimumFeatureSubsetSize){ break; }
    }
    
//    DEBUGMSG_VAR(isSelected);
    for (unsigned i = 0; i < numberOfTotalFeatures; i++) {
      if (isSelected[i] == true){ featureSubsetIndices.push_back(i); }
    }
//    DEBUGMSG_VAR(numberOfTotalFeatures);
//    DEBUGMSG_VAR(featureSubsetIndices);
    return featureSubsetIndices;
  }
  
    // TODO: we need to modify this later for randomized random forest algo's implementation
  int getOptimumFeatureSubsetSize(){ return (int) ceil(log(numberOfTotalFeatures) / log(2)); }
  
  void createBootstrappedSamples(){
    
      // randomly create indices list
    vector<bool> isInTrainingSamples;
    
    for (unsigned i = 0; i < bootstrappedSamplesSize; i++) {
      isInTrainingSamples.push_back(false);
    }
    
    bootstrappedTrainingSamples.clear();
    for (unsigned i = 0; i < bootstrappedSamplesSize; i++) {
      int randomIndex = (int)(((double)(rand()) / (double)(RAND_MAX)) * bootstrappedSamplesSize);
//      DEBUGMSG_VAR(randomIndex);
      isInTrainingSamples[randomIndex] = true;
      bootstrappedTrainingSamples.push_back(baseSamples[randomIndex]);
    }
    
    bootStrammpedTestSamples.clear();
    for (unsigned i = 0; i < bootstrappedSamplesSize; i++) {
      if(isInTrainingSamples[i] == false){
        bootStrammpedTestSamples.push_back(baseSamples[i]);
      }
    }
    
//    DEBUGMSG_VAR(bootStrammpedTestSamples.size());
    
  }
  
    // function that is useful for debugging
  void printSamples(vector<TrainingSet> samples){
    for (unsigned i = 0; i < samples.size(); i++) {
      vector<int> otuCounts = samples[i].getOtuCounts();
      int outputClassId = samples[i].getOutputClassId();
      string outputClass = samples[i].getOutputClass();
      for (unsigned j = 0; j < otuCounts.size(); j++) {
        cout << otuCounts[j] << " ";
      }
      cout << outputClass << " " << outputClassId << endl;
    }
  }
  
private:
  vector<TrainingSet> baseSamples;
  vector<TrainingSet> bootstrappedTrainingSamples;
  vector<TrainingSet> bootStrammpedTestSamples;
  
  int bootstrappedSamplesSize;
  int numberOfTotalFeatures;
};


#endif
