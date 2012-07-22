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
  splitFeatureIndex(0),
  splitFeatureValue(0),
  splitFeatureEntropy(0),
  ownEntropy(0),
  leftChildNode(NULL),
  rightChildNode(NULL),
  parentNode(NULL){
//    self.bootstrappedFeatureVectors = [list(x) for x in zip(*self.bootstrappedTrainingSamples)]
//    self.bootstrappedOutputVector = [bootstrappedTrainingSamples[x][self.numFeatures] for x in range(0, self.numSamples)]
//    
//# call some helper functions
//    self.createLocalDiscardedFeatureList()
//    
//    self.updateNodeEntropy()

  }
  
  ~TreeNode(){
  }
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
  int splitFeatureEntropy;
  int ownEntropy;
  
  vector<vector<int> > bootstrappedFeatureVectors;
  vector<int> bootstrappedOutputVector;
  vector<int> featureSubsetIndices;
  
  TreeNode* leftChildNode;
  TreeNode* rightChildNode;
  TreeNode* parentNode;
  
  vector<int> localDiscardedFeatureIndices;
  
  void createLocalDiscardedFeatureList(){
  }
  
  void updateNodeEntropy(){
  }
};

#endif
