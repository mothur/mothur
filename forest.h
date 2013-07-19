//
//  forest.h
//  Mothur
//
//  Created by Kathryn Iverson on 10/26/12. Modified abstractrandomforest
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__forest__
#define __Mothur__forest__

#include <iostream>
#include "mothurout.h"
#include "macros.h"
#include "decisiontree.hpp"
#include "abstractdecisiontree.hpp"
/***********************************************************************/
//this is a re-implementation of the abstractrandomforest class

class Forest{
public:
    // intialization with vectors
    Forest(const std::vector < std::vector<int> > dataSet,
           const int numDecisionTrees,
           const string treeSplitCriterion,
           const bool doPruning,
           const float pruneAggressiveness,
           const bool discardHighErrorTrees,
           const float highErrorTreeDiscardThreshold,
           const string optimumFeatureSubsetSelectionCriteria,
           const float featureStandardDeviationThreshold);
    virtual ~Forest(){ }
    virtual int populateDecisionTrees() = 0;
    virtual int calcForrestErrorRate() = 0;
    virtual int calcForrestVariableImportance(string) = 0;
    virtual int updateGlobalOutOfBagEstimates(DecisionTree* decisionTree) = 0;
    
    /***********************************************************************/
    
protected:
    
    // TODO: create a better way of discarding feature
    // currently we just set FEATURE_DISCARD_SD_THRESHOLD to 0 to solved this
    // it can be tuned for better selection
    // also, there might be other factors like Mean or other stuffs
    // same would apply for createLocalDiscardedFeatureList in the TreeNode class
    
    // TODO: Another idea is getting an aggregated discarded feature indices after the run, from combining
    // the local discarded feature indices
    // this would penalize a feature, even if in global space the feature looks quite good
    // the penalization would be averaged, so this woould unlikely to create a local optmina
    
    vector<int> getGlobalDiscardedFeatureIndices();
    void calculateFScore();
    
    int numDecisionTrees;
    int numSamples;
    int numFeatures;
    vector< vector<int> > dataSet;
    vector<int> globalDiscardedFeatureIndices;
    vector<double> globalVariableImportanceList;
    
    vector< pair<int, double> > featureRanksByFScore;
    
    string treeSplitCriterion;
  
    bool doPruning;
    float pruneAggressiveness;
    bool discardHighErrorTrees;
    float highErrorTreeDiscardThreshold;
    string optimumFeatureSubsetSelectionCriteria;
    float featureStandardDeviationThreshold;
  
    // This is a map of each feature to outcome count of each classes
    // e.g. 1 => [2 7] means feature 1 has 2 outcome of 0 and 7 outcome of 1
    map<int, vector<int> > globalOutOfBagEstimates;
    
    // TODO: fix this, do we use pointers?
    vector<AbstractDecisionTree*> decisionTrees;
    
    // predictedClasses[i] denotes the class predicted by the algorithm
    // for i'th training sample
    vector<int> predictedClasses;

    MothurOut* m;
    
private:
    
};

#endif /* defined(__Mothur__forest__) */
