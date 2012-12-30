//
//  regularizedrandomforest.h
//  Mothur
//
//  Created by Kathryn Iverson on 11/16/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#ifndef __Mothur__regularizedrandomforest__
#define __Mothur__regularizedrandomforest__

#include "forest.h"
#include "decisiontree.hpp"

class RegularizedRandomForest: public Forest {
public:
    //
    RegularizedRandomForest(const vector <vector<int> > dataSet,
                            const int numDecisionTrees,
                            const string treeSplitCriterion,
                            const bool doPruning,
                            const float pruneAggressiveness,
                            const bool discardHighErrorTrees,
                            const float highErrorTreeDiscardThreshold,
                            const string optimumFeatureSubsetSelectionCriteria,
                            const float featureStandardDeviationThreshold);
    
    int calcForrestErrorRate();
    int calcForrestVariableImportance(string);
    int populateDecisionTrees();
    int updateGlobalOutOfBagEstimates(DecisionTree* decisionTree);
    
private:
    //
    MothurOut* m;
};

#endif /* defined(__Mothur__regularizedrandomforest__) */
