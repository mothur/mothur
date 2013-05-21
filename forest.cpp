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
               const string treeSplitCriterion = "gainratio",
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
        globalVariableImportanceList(numFeatures, 0),
        treeSplitCriterion(treeSplitCriterion),
        doPruning(doPruning),
        pruneAggressiveness(pruneAggressiveness),
        discardHighErrorTrees(discardHighErrorTrees),
        highErrorTreeDiscardThreshold(highErrorTreeDiscardThreshold),
        optimumFeatureSubsetSelectionCriteria(optimumFeatureSubsetSelectionCriteria),
        featureStandardDeviationThreshold(featureStandardDeviationThreshold)
        {
        
    m = MothurOut::getInstance();
    globalDiscardedFeatureIndices = getGlobalDiscardedFeatureIndices();
//    calculateFScore();
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
            if (standardDeviation <= featureStandardDeviationThreshold){ globalDiscardedFeatureIndices.push_back(i); }
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

void Forest::calculateFScore() {
    vector< vector<int> > featureVectors(numFeatures, vector<int>(numSamples, 0) );
    vector<int> outputVector(numSamples, 0);
    vector<int> outputClasses;
    int numOutputClasses = 0;

        // calculate feature vectors and output vector
    for (int i = 0; i < numSamples; i++) {
        for (int j = 0; j < numFeatures; j++) { featureVectors[j][i] = dataSet[i][j]; }
        int outcome = dataSet[i][numFeatures];
        outputVector[i] = outcome;
        vector<int>::iterator it = find(outputClasses.begin(), outputClasses.end(), outcome);
        if (it == outputClasses.end()){       // find() will return classes.end() if the element is not found
            outputClasses.push_back(outcome);
            numOutputClasses++;
        }
    }    
    
    for (int i = 0; i < numFeatures; i++) {
        double fisherNumerator = 0;
        double fisherDenominator = 0;
        double fisherScore = 0;
        
        long featureSum = accumulate(featureVectors[i].begin(), featureVectors[i].end(), 0);
        double featureMean = (double) featureSum / (double) numSamples;
        
        vector<int> outputClassCounts(numOutputClasses, 0);
        vector<long> outputClassBasedSums(numOutputClasses, 0);
        vector<double> outputClassBasedMeans(numOutputClasses, 0);
        
        for (int j = 0; j < numSamples; j++) {
            outputClassCounts[outputVector[j]]++;
            outputClassBasedSums[outputVector[j]] += featureVectors[i][j];
        }
        
        for (int k = 0; k < numOutputClasses; k++) {
            outputClassBasedMeans[k] = (double) outputClassBasedSums[k] / (double)outputClassCounts[k];
            fisherNumerator += pow(outputClassBasedMeans[k] - featureMean, 2);
        }
        
        for (int j = 0; j < numSamples; j++) {
                // TODO: choose between weighted and un-weighted denominator calculation
                // which one is the correct and/or perfect?
            
                // This is the unweighted calculation
//            fisherDenominator += pow(featureVectors[i][j] - outputClassBasedMeans[outputVector[j]], 2);
            
                // this is the weighted calculation
            fisherDenominator += (1 / ((double)outputClassCounts[outputVector[j]] - 1)) * pow(featureVectors[i][j] - outputClassBasedMeans[outputVector[j]], 2);
        }
        
        if (fisherDenominator == 0) { fisherDenominator = numeric_limits<double>::min(); }
        
        fisherScore = fisherNumerator / fisherDenominator;
        pair<int, double> fisherScoreRank(i, fisherScore);
        featureRanksByFScore.push_back(fisherScoreRank);
    }
    
    VariableRankDescendingSorterDouble sorter;
    sort(featureRanksByFScore.begin(), featureRanksByFScore.end(), sorter);
        
    // scale the scores, the top feature will have it's value scaled to 1
//    for (vector<pair<int, double> >::iterator it = featureRanksByFScore.begin(); it != featureRanksByFScore.end(); it++) {
//        it->second /= featureRanksByFScore.begin()->second;
//        if (it->second < 0) {
//            featureRanksByFScore.erase(it, featureRanksByFScore.end());
//            break;
//        }
//    }
    
//    for (int i = 0; i < featureRanksByFScore.size(); i++) {
//        cout << featureRanksByFScore[i].first << "\t" << featureRanksByFScore[i].second << endl;
//    }
//    cout << endl;
//    exit(1);
}

