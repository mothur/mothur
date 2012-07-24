//
//  treenode.hpp
//  rrf-fs-prototype
//
//  Created by Abu Zaher Faridee on 5/29/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef rrf_fs_prototype_treenode_hpp
#define rrf_fs_prototype_treenode_hpp

class TreeNode{
public:
  TreeNode(vector< vector<int> > bootstrappedTrainingSamples,
           vector<int> globalDiscardedFeatureIndices,
           int numFeatures,
           int numSamples,
           int numOutputClasses,
           int generation)
  : bootstrappedTrainingSamples(bootstrappedTrainingSamples),
  globalDiscardedFeatureIndices(globalDiscardedFeatureIndices),
  numFeatures(numFeatures),
  numSamples(numSamples),
  numOutputClasses(numOutputClasses),
  generation(generation),
  isLeaf(false),
  outputClass(-1),
  splitFeatureIndex(-1),
  splitFeatureValue(-1),
  splitFeatureEntropy(-1.0),
  ownEntropy(-1.0),
  bootstrappedFeatureVectors(numFeatures, vector<int>(numSamples, 0)),
  bootstrappedOutputVector(numSamples, 0),
  leftChildNode(NULL),
  rightChildNode(NULL),
  parentNode(NULL){

    for (int i = 0; i < numSamples; i++) {    // just doing a simple transpose of the matrix
      for (int j = 0; j < numFeatures; j++) { bootstrappedFeatureVectors[j][i] = bootstrappedTrainingSamples[i][j]; }
    }
    
    for (int i = 0; i < numSamples; i++) { bootstrappedOutputVector[i] = bootstrappedTrainingSamples[i][numFeatures]; }
    
    createLocalDiscardedFeatureList();
    updateNodeEntropy();
  }
  
  ~TreeNode(){
  }
  
    // getters
    // we need to return const reference so that we have the actual value and note a copy, 
    // plus we do not modify the value as well
  const int& getSplitFeatureIndex() { return splitFeatureIndex; }
    // TODO: check if this works properly or returs a shallow copy of the data
  const vector< vector<int> >& getBootstrappedTrainingSamples() { return bootstrappedTrainingSamples; }
  const int& getSplitFeatureValue() { return splitFeatureValue; }
  const int& getGeneration() { return generation; }
  const bool& checkIsLeaf() { return isLeaf; }
    // TODO: fix this const pointer dillema
    // we do not want to modify the data pointer by getLeftChildNode
  TreeNode* getLeftChildNode() { return leftChildNode; }
  TreeNode* getRightChildNode() { return rightChildNode; }
  const int& getOutputClass() { return outputClass; }
  const int& getNumSamples() { return numSamples; }
  const int& getNumFeatures() { return numFeatures; }
  const vector<int>& getLocalDiscardedFeatureIndices() { return localDiscardedFeatureIndices; }
  const vector< vector<int> >& getBootstrappedFeatureVectors() { return bootstrappedFeatureVectors; }
  const vector<int>& getBootstrappedOutputVector() { return bootstrappedOutputVector; }
  const vector<int>& getFeatureSubsetIndices() { return featureSubsetIndices; }
  const double& getOwnEntropy() { return ownEntropy; }
  
    // setters
  void setIsLeaf(bool isLeaf) { this->isLeaf = isLeaf; }
  void setOutputClass(int outputClass) { this->outputClass = outputClass; }
  void setFeatureSubsetIndices(vector<int> featureSubsetIndices) { this->featureSubsetIndices = featureSubsetIndices; }
  void setLeftChildNode(TreeNode* leftChildNode) { this->leftChildNode = leftChildNode; }
  void setRightChildNode(TreeNode* rightChildNode) { this->rightChildNode = rightChildNode; }
  void setParentNode(TreeNode* parentNode) { this->parentNode = parentNode; }
  void setSplitFeatureIndex(int splitFeatureIndex) { this->splitFeatureIndex = splitFeatureIndex; }
  void setSplitFeatureValue(int splitFeatureValue) { this->splitFeatureValue = splitFeatureValue; }
  void setSplitFeatureEntropy(double splitFeatureEntropy) { this->splitFeatureEntropy = splitFeatureEntropy; }
  
protected:
private:
  vector<vector<int> > bootstrappedTrainingSamples;
  vector<int> globalDiscardedFeatureIndices;
  int numFeatures;
  int numSamples;
  int numOutputClasses;
  int generation;
  bool isLeaf;
  int outputClass;
  
  int splitFeatureIndex;
  int splitFeatureValue;
  double splitFeatureEntropy;
  double ownEntropy;
  
  vector<vector<int> > bootstrappedFeatureVectors;
  vector<int> bootstrappedOutputVector;
  vector<int> featureSubsetIndices;
  
  TreeNode* leftChildNode;
  TreeNode* rightChildNode;
  TreeNode* parentNode;
  
  vector<int> localDiscardedFeatureIndices;
  
  void createLocalDiscardedFeatureList(){
    for (int i = 0; i < numFeatures; i++) {
      vector<int>::iterator it = find(globalDiscardedFeatureIndices.begin(), globalDiscardedFeatureIndices.end(), i);
      if (it == globalDiscardedFeatureIndices.end()){                           // NOT FOUND
        double standardDeviation = getStandardDeviation(bootstrappedFeatureVectors[i]);  
        if (standardDeviation <= 0){ localDiscardedFeatureIndices.push_back(i); }
      }
    }
  }
  
  void updateNodeEntropy(){    
    vector<int> classCounts(numOutputClasses, 0);
    for (int i = 0; i < bootstrappedOutputVector.size(); i++) { classCounts[bootstrappedOutputVector[i]]++; }
    int totalClassCounts = accumulate(classCounts.begin(), classCounts.end(), 0);
    double nodeEntropy = 0.0;
    for (int i = 0; i < classCounts.size(); i++) {
      if (classCounts[i] == 0) continue;
      double probability = (double)classCounts[i] / (double)totalClassCounts;
      nodeEntropy += -(probability * log2(probability));
    }
    ownEntropy = nodeEntropy;
  }
  
};

#endif
