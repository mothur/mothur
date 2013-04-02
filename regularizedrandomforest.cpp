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
                                                 const string treeSplitCriterion = "gainratio")
                      : Forest(dataSet,
                               numDecisionTrees,
                               treeSplitCriterion,
                               false, 0.9, true, 0.4, "log2", 0.0) {
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
