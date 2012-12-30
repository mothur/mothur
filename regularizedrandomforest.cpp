//
//  regularizedrandomforest.cpp
//  Mothur
//
//  Created by Kathryn Iverson on 11/16/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "regularizedrandomforest.h"

RegularizedRandomForest::RegularizedRandomForest(const vector <vector<int> > dataSet,
                           const int numDecisionTrees,
                           const string treeSplitCriterion = "informationGain",
                           const bool doPruning = false,
                           const float pruneAggressiveness = 0.9,
                           const bool discardHighErrorTrees = true,
                           const float highErrorTreeDiscardThreshold = 0.4,
                           const string optimumFeatureSubsetSelectionCriteria = "log2",
                           const float featureStandardDeviationThreshold = 0.0)
: Forest(dataSet, numDecisionTrees, treeSplitCriterion, doPruning, pruneAggressiveness, discardHighErrorTrees, highErrorTreeDiscardThreshold, optimumFeatureSubsetSelectionCriteria, featureStandardDeviationThreshold) {
    m = MothurOut::getInstance();
}

int RegularizedRandomForest::calcForrestErrorRate() {
    //
    try {
        return 0;
    }
    catch(exception& e) {
		m->errorOut(e, "RegularizedRandomForest", "calcForrestErrorRate");
		exit(1);
	}
}

int RegularizedRandomForest::calcForrestVariableImportance(string filename) {
    //
    try {
        return 0;
    }
    catch(exception& e) {
		m->errorOut(e, "RegularizedRandomForest", "calcForrestVariableImportance");
		exit(1);
	}
}

int RegularizedRandomForest::populateDecisionTrees() {
    //
    try {
        return 0;
    }
    catch(exception& e) {
		m->errorOut(e, "RegularizedRandomForest", "populateDecisionTrees");
		exit(1);
	}
}

int RegularizedRandomForest::updateGlobalOutOfBagEstimates(DecisionTree *decisionTree) {
    try {
        return 0;
    }
    catch(exception& e) {
		m->errorOut(e, "RegularizedRandomForest", "updateGlobalOutOfBagEstimates");
		exit(1);
	}
}