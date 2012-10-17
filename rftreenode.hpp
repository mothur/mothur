//
//  rftreenode.hpp
//  rrf-fs-prototype
//
//  Created by Abu Zaher Faridee on 5/29/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#ifndef rrf_fs_prototype_treenode_hpp
#define rrf_fs_prototype_treenode_hpp

#include "mothurout.h"
#include "macros.h"

class RFTreeNode{
    
public:
    
    RFTreeNode(vector< vector<int> > bootstrappedTrainingSamples, vector<int> globalDiscardedFeatureIndices, int numFeatures, int numSamples, int numOutputClasses, int generation);
    
    virtual ~RFTreeNode(){}
    
    // getters
    // we need to return const reference so that we have the actual value and not a copy, 
    // plus we do not modify the value as well
    const int getSplitFeatureIndex() { return splitFeatureIndex; }
    // TODO: check if this works properly or returs a shallow copy of the data
    const vector< vector<int> >& getBootstrappedTrainingSamples() { return bootstrappedTrainingSamples; }
    const int getSplitFeatureValue() { return splitFeatureValue; }
    const int getGeneration() { return generation; }
    const bool checkIsLeaf() { return isLeaf; }
    // TODO: fix this const pointer dillema
    // we do not want to modify the data pointer by getLeftChildNode
    RFTreeNode* getLeftChildNode() { return leftChildNode; }
    RFTreeNode* getRightChildNode() { return rightChildNode; }
    const int getOutputClass() { return outputClass; }
    const int getNumSamples() { return numSamples; }
    const int getNumFeatures() { return numFeatures; }
    const vector<int>& getLocalDiscardedFeatureIndices() { return localDiscardedFeatureIndices; }
    const vector< vector<int> >& getBootstrappedFeatureVectors() { return bootstrappedFeatureVectors; }
    const vector<int>& getBootstrappedOutputVector() { return bootstrappedOutputVector; }
    const vector<int>& getFeatureSubsetIndices() { return featureSubsetIndices; }
    const double getOwnEntropy() { return ownEntropy; }
    
    // setters
    void setIsLeaf(bool isLeaf) { this->isLeaf = isLeaf; }
    void setOutputClass(int outputClass) { this->outputClass = outputClass; }
    void setFeatureSubsetIndices(vector<int> featureSubsetIndices) { this->featureSubsetIndices = featureSubsetIndices; }
    void setLeftChildNode(RFTreeNode* leftChildNode) { this->leftChildNode = leftChildNode; }
    void setRightChildNode(RFTreeNode* rightChildNode) { this->rightChildNode = rightChildNode; }
    void setParentNode(RFTreeNode* parentNode) { this->parentNode = parentNode; }
    void setSplitFeatureIndex(int splitFeatureIndex) { this->splitFeatureIndex = splitFeatureIndex; }
    void setSplitFeatureValue(int splitFeatureValue) { this->splitFeatureValue = splitFeatureValue; }
    void setSplitFeatureEntropy(double splitFeatureEntropy) { this->splitFeatureEntropy = splitFeatureEntropy; }
    
    // TODO: need to remove this mechanism of friend class
    //NOTE: friend classes can be useful for testing purposes, but I would avoid using them otherwise.
    friend class DecisionTree;
    friend class AbstractDecisionTree;
    
private:
    vector<vector<int> > bootstrappedTrainingSamples;
    vector<int> globalDiscardedFeatureIndices;
    vector<int> localDiscardedFeatureIndices;
    vector<vector<int> > bootstrappedFeatureVectors;
    vector<int> bootstrappedOutputVector;
    vector<int> featureSubsetIndices;

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
    
    RFTreeNode* leftChildNode;
    RFTreeNode* rightChildNode;
    RFTreeNode* parentNode;
    
    MothurOut* m;
    
    int createLocalDiscardedFeatureList();
    int updateNodeEntropy();
    
};

#endif
