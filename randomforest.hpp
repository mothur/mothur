//
//  randomforest.hpp
//  rrf-fs-prototype
//
//  Created by Abu Zaher Faridee on 7/20/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#ifndef RF_RANDOMFOREST_HPP
#define RF_RANDOMFOREST_HPP

#include "macros.h"
#include "forest.h"
#include "decisiontree.hpp"

class RandomForest: public Forest {
    
public:
    
    // DONE
    RandomForest(const vector <vector<int> > dataSet,
                 const int numDecisionTrees,
                 const string treeSplitCriterion,
                 const bool doPruning,
                 const float pruneAggressiveness,
                 const bool discardHighErrorTrees,
                 const float highErrorTreeDiscardThreshold,
                 const string optimumFeatureSubsetSelectionCriteria,
                 const float featureStandardDeviationThreshold);
    
    
    //NOTE:: if you are going to dynamically cast, aren't you undoing the advantage of abstraction. Why abstract at all?
    //could cause maintenance issues later if other types of Abstract decison trees are created that cannot be cast as a decision tree.
//    virtual ~RandomForest() {
//        for (vector<AbstractDecisionTree*>::iterator it = decisionTrees.begin(); it != decisionTrees.end(); it++) {
//            // we know that this is decision tree, so we can do a dynamic_case<DecisionTree*> here
//            DecisionTree* decisionTree = dynamic_cast<DecisionTree*>(*it);
//            // calling the destructor by deleting
//            delete decisionTree;
//        }
//    }
    
    int calcForrestErrorRate();
    int calcForrestVariableImportance(string);
    int populateDecisionTrees();
    int updateGlobalOutOfBagEstimates(DecisionTree* decisionTree);
    
private:
    MothurOut* m;
    
};

#endif
