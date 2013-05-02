//
//  abstractdecisiontree.hpp
//  rrf-fs-prototype
//
//  Created by Abu Zaher Faridee on 7/22/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#ifndef RF_ABSTRACTDECISIONTREE_HPP
#define RF_ABSTRACTDECISIONTREE_HPP

#include "mothurout.h"
#include "macros.h"
#include "rftreenode.hpp"

#define DEBUG_MODE

/**************************************************************************************************/

struct IntPairVectorSorter{
    bool operator() (const pair<int, int>& firstPair, const pair<int, int>& secondPair) {
        return firstPair.first < secondPair.first;
    }
};

/**************************************************************************************************/

class AbstractDecisionTree{
  
public:
  
    AbstractDecisionTree(vector<vector<int> >& baseDataSet,
                           vector<int> globalDiscardedFeatureIndices, 
                           OptimumFeatureSubsetSelector optimumFeatureSubsetSelector, 
                           string treeSplitCriterion);    
    virtual ~AbstractDecisionTree(){}
    
  
protected:
  
    virtual int createBootStrappedSamples();
    virtual int getMinEntropyOfFeature(vector<int> featureVector, vector<int> outputVector, double& minEntropy, int& featureSplitValue, double& intrinsicValue);
        virtual int getBestSplitAndMinEntropy(vector< pair<int, int> > featureOutputPairs, vector<int> splitPoints, double& minEntropy, int& minEntropyIndex, double& relatedIntrinsicValue);
    virtual double calcIntrinsicValue(int numLessThanValueAtSplitPoint, int numGreaterThanValueAtSplitPoint, int numSamples);
    virtual double calcSplitEntropy(vector< pair<int, int> > featureOutputPairs, int splitIndex, int numOutputClasses, bool);

    virtual int getSplitPopulation(RFTreeNode* node, vector< vector<int> >& leftChildSamples, vector< vector<int> >& rightChildSamples);
    virtual bool checkIfAlreadyClassified(RFTreeNode* treeNode, int& outputClass);

    vector< vector<int> >& baseDataSet;
    int numSamples;
    int numFeatures;
    int numOutputClasses;
    vector<int> outputClasses;
    
    vector< vector<int> > bootstrappedTrainingSamples;
    vector<int> bootstrappedTrainingSampleIndices;
    vector< vector<int> > bootstrappedTestSamples;
    vector<int> bootstrappedTestSampleIndices;
    
    vector<vector<int> > testSampleFeatureVectors;
    
    RFTreeNode* rootNode;
    int nodeIdCount;
    map<int, int> nodeMisclassificationCounts;
    vector<int> globalDiscardedFeatureIndices;
    int optimumFeatureSubsetSize;
    string treeSplitCriterion;
    MothurOut* m;
  
private:
    
  
};
/**************************************************************************************************/

#endif
