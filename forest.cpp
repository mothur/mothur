//
//  forest.cpp
//  Mothur
//
//  Created by Kathryn Iverson on 10/26/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "forest.h"

/***********************************************************************/
Forest::Forest(const std::vector < std::vector<int> > dataSet,
                                           const int numDecisionTrees,
                                           const string treeSplitCriterion = "informationGain")
: dataSet(dataSet),
numDecisionTrees(numDecisionTrees),
numSamples((int)dataSet.size()),
numFeatures((int)(dataSet[0].size() - 1)),
globalVariableImportanceList(numFeatures, 0),
treeSplitCriterion(treeSplitCriterion) {
    m = MothurOut::getInstance();
    globalDiscardedFeatureIndices = getGlobalDiscardedFeatureIndices();
    // TODO: double check if the implemenatation of 'globalOutOfBagEstimates' is correct
}

/***********************************************************************/
Forest::Forest(const std::vector < std::vector<int> > dataSet,
                                           const int numDecisionTrees,
                                           const string treeSplitCriterion = "informationGain",
                                           const bool doPruning = false,
                                           const float pruneAggressiveness = 0.9,
                                           const bool discardHighErrorTrees = true,
                                           const float highErrorTreeDiscardThreshold = 0.4,
                                           const string optimumFeatureSubsetSelectionCriteria = "log2",
                                           const float featureStandardDeviationThreshold = 0.0)
: dataSet(dataSet),
numDecisionTrees(numDecisionTrees),
numSamples((int)dataSet.size()),
numFeatures((int)(dataSet[0].size() - 1)),
globalDiscardedFeatureIndices(getGlobalDiscardedFeatureIndices()),
globalVariableImportanceList(numFeatures, 0),
treeSplitCriterion(treeSplitCriterion),
doPruning(doPruning),
pruneAggressiveness(pruneAggressiveness),
discardHighErrorTrees(discardHighErrorTrees),
highErrorTreeDiscardThreshold(highErrorTreeDiscardThreshold),
optimumFeatureSubsetSelectionCriteria(optimumFeatureSubsetSelectionCriteria),
featureStandardDeviationThreshold(featureStandardDeviationThreshold) {
    m = MothurOut::getInstance();
        // TODO: double check if the implemenatation of 'globalOutOfBagEstimates' is correct
}

/***********************************************************************/

vector<int> Forest::getGlobalDiscardedFeatureIndices() {
    try {
        //vector<int> globalDiscardedFeatureIndices;
        //globalDiscardedFeatureIndices.push_back(1);
        
        // calculate feature vectors
        vector< vector<int> > featureVectors(numFeatures, vector<int>(numSamples, 0) );
        for (int i = 0; i < numSamples; i++) {
            if (m->control_pressed) { return globalDiscardedFeatureIndices; }
            for (int j = 0; j < numFeatures; j++) { featureVectors[j][i] = dataSet[i][j]; }
        }
        
        for (int i = 0; i < featureVectors.size(); i++) {
            if (m->control_pressed) { return globalDiscardedFeatureIndices; }
            double standardDeviation = m->getStandardDeviation(featureVectors[i]);
            if (standardDeviation <= 0){ globalDiscardedFeatureIndices.push_back(i); }
        }
        
        if (m->debug) {
            m->mothurOut("number of global discarded features:  " + toString(globalDiscardedFeatureIndices.size())+ "\n");
            m->mothurOut("total features: " + toString(featureVectors.size())+ "\n");
        }
        
        return globalDiscardedFeatureIndices;
    }
	catch(exception& e) {
		m->errorOut(e, "Forest", "getGlobalDiscardedFeatureIndices");
		exit(1);
	}
}

/***********************************************************************/

