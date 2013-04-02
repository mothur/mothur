//
//  randomforest.cpp
//  Mothur
//
//  Created by Sarah Westcott on 10/2/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "randomforest.hpp" 

/***********************************************************************/

RandomForest::RandomForest(const vector <vector<int> > dataSet,
                           const int numDecisionTrees,
                           const string treeSplitCriterion = "gainratio",
                           const bool doPruning = false,
                           const float pruneAggressiveness = 0.9,
                           const bool discardHighErrorTrees = true,
                           const float highErrorTreeDiscardThreshold = 0.4,
                           const string optimumFeatureSubsetSelectionCriteria = "log2",
                           const float featureStandardDeviationThreshold = 0.0)
            : Forest(dataSet, numDecisionTrees, treeSplitCriterion, doPruning, pruneAggressiveness, discardHighErrorTrees, highErrorTreeDiscardThreshold, optimumFeatureSubsetSelectionCriteria, featureStandardDeviationThreshold) {
    m = MothurOut::getInstance();
}

/***********************************************************************/
// DONE
int RandomForest::calcForrestErrorRate() {
    try {
        int numCorrect = 0;
        for (map<int, vector<int> >::iterator it = globalOutOfBagEstimates.begin(); it != globalOutOfBagEstimates.end(); it++) {
            
            if (m->control_pressed) { return 0; }
            
            int indexOfSample = it->first;
            vector<int> predictedOutComes = it->second;
            vector<int>::iterator maxPredictedOutComeIterator = max_element(predictedOutComes.begin(), predictedOutComes.end());
            int majorityVotedOutcome = (int)(maxPredictedOutComeIterator - predictedOutComes.begin());
            int realOutcome = dataSet[indexOfSample][numFeatures];
            
            if (majorityVotedOutcome == realOutcome) { numCorrect++; }
        }
        
        // TODO: save or return forrestErrorRate for future use;
        double forrestErrorRate = 1 - ((double)numCorrect / (double)globalOutOfBagEstimates.size());
        
        m->mothurOut("numCorrect = " + toString(numCorrect)+ "\n");
        m->mothurOut("forrestErrorRate = " + toString(forrestErrorRate)+ "\n");
    
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "RandomForest", "calcForrestErrorRate");
		exit(1);
	} 
}

/***********************************************************************/
// DONE
int RandomForest::calcForrestVariableImportance(string filename) {
    try {
    
        // follow the link: http://en.wikipedia.org/wiki/Dynamic_cast
        //if you are going to dynamically cast, aren't you undoing the advantage of abstraction. Why abstract at all?
        //could cause maintenance issues later if other types of Abstract decison trees are created that cannot be cast as a decision tree.
        for (int i = 0; i < decisionTrees.size(); i++) {
            if (m->control_pressed) { return 0; }
            
            DecisionTree* decisionTree = dynamic_cast<DecisionTree*>(decisionTrees[i]);
            
            for (int j = 0; j < numFeatures; j++) {
                globalVariableImportanceList[j] += (double)decisionTree->variableImportanceList[j];
            }
        }
        
        for (int i = 0;  i < numFeatures; i++) {
            cout << "[" << i << ',' << globalVariableImportanceList[i] << "], ";
            globalVariableImportanceList[i] /= (double)numDecisionTrees;
        }
        
        vector< vector<double> > globalVariableRanks;
        for (int i = 0; i < globalVariableImportanceList.size(); i++) {
            if (globalVariableImportanceList[i] > 0) {
                vector<double> globalVariableRank(2, 0);
                globalVariableRank[0] = i; globalVariableRank[1] = globalVariableImportanceList[i];
                globalVariableRanks.push_back(globalVariableRank);
            }
        }
        
        VariableRankDescendingSorterDouble variableRankDescendingSorter;
        sort(globalVariableRanks.begin(), globalVariableRanks.end(), variableRankDescendingSorter);
        
        ofstream out;
        m->openOutputFile(filename, out);
        out <<"OTU\tRank\n";
        for (int i = 0; i < globalVariableRanks.size(); i++) {
            out << m->currentBinLabels[(int)globalVariableRanks[i][0]] << '\t' << globalVariableImportanceList[globalVariableRanks[i][0]] << endl;
        }
        out.close();
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "RandomForest", "calcForrestVariableImportance");
		exit(1);
	}  
}
/***********************************************************************/
int RandomForest::populateDecisionTrees() {
    try {
        
        vector<double> errorRateImprovements;
        
        for (int i = 0; i < numDecisionTrees; i++) {
          
            if (m->control_pressed) { return 0; }
            if (((i+1) % 10) == 0) {  m->mothurOut("Creating " + toString(i+1) + " (th) Decision tree\n");  }
          
            // TODO: need to first fix if we are going to use pointer based system or anything else
            DecisionTree* decisionTree = new DecisionTree(dataSet, globalDiscardedFeatureIndices, OptimumFeatureSubsetSelector(optimumFeatureSubsetSelectionCriteria), treeSplitCriterion, featureStandardDeviationThreshold);
          
            if (m->debug && doPruning) {
                m->mothurOut("Before pruning\n");
                decisionTree->printTree(decisionTree->rootNode, "ROOT");
            }
            
            int numCorrect;
            double treeErrorRate;
            
            decisionTree->calcTreeErrorRate(numCorrect, treeErrorRate);
            double prePrunedErrorRate = treeErrorRate;
            
            if (m->debug) {
                m->mothurOut("treeErrorRate: " + toString(treeErrorRate) + " numCorrect: " + toString(numCorrect) + "\n");
            }
            
            if (doPruning) {
                decisionTree->pruneTree(pruneAggressiveness);
                if (m->debug) {
                    m->mothurOut("After pruning\n");
                    decisionTree->printTree(decisionTree->rootNode, "ROOT");
                }
                decisionTree->calcTreeErrorRate(numCorrect, treeErrorRate);
            }
            double postPrunedErrorRate = treeErrorRate;
            
          
            decisionTree->calcTreeVariableImportanceAndError(numCorrect, treeErrorRate);
            if (m->debug && doPruning) {
                m->mothurOut("treeErrorRate: " + toString(treeErrorRate) + " numCorrect: " + toString(numCorrect) + "\n");
            }
            
            double errorRateImprovement = (prePrunedErrorRate - postPrunedErrorRate) / prePrunedErrorRate;
            m->mothurOut("errorRateImprovement: " + toString(errorRateImprovement) + "\n");
                        
            if (discardHighErrorTrees) {
                if (treeErrorRate < highErrorTreeDiscardThreshold) {
                    updateGlobalOutOfBagEstimates(decisionTree);
                    decisionTree->purgeDataSetsFromTree();
                    decisionTrees.push_back(decisionTree);
                    if (doPruning) {
                        errorRateImprovements.push_back(errorRateImprovement);
                    }
                } else {
                    delete decisionTree;
                }
            } else {
                updateGlobalOutOfBagEstimates(decisionTree);
                decisionTree->purgeDataSetsFromTree();
                decisionTrees.push_back(decisionTree);
                if (doPruning) {
                    errorRateImprovements.push_back(errorRateImprovement);
                }
            }          
        }
        
        double avgErrorRateImprovement = -1.0;
        if (errorRateImprovements.size() > 0) {
            avgErrorRateImprovement = accumulate(errorRateImprovements.begin(), errorRateImprovements.end(), 0.0);
            cout << "Total " << avgErrorRateImprovement << " size " << errorRateImprovements.size() << endl;
            avgErrorRateImprovement /= errorRateImprovements.size();
        }
        
        if (m->debug) {
            m->mothurOut("avgErrorRateImprovement:" + toString(avgErrorRateImprovement) + "\n");
            // m->mothurOut("globalOutOfBagEstimates = " + toStringVectorMap(globalOutOfBagEstimates)+ "\n");
        }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "RandomForest", "populateDecisionTrees");
        exit(1);
    }  
}
/***********************************************************************/
// TODO: need to finalize bettween reference and pointer for DecisionTree [partially solved]
// DONE: make this pure virtual in superclass
// DONE
int RandomForest::updateGlobalOutOfBagEstimates(DecisionTree* decisionTree) {
    try {
        for (map<int, int>::iterator it = decisionTree->outOfBagEstimates.begin(); it != decisionTree->outOfBagEstimates.end(); it++) {
            
            if (m->control_pressed) { return 0; }
            
            int indexOfSample = it->first;
            int predictedOutcomeOfSample = it->second;
            
            if (globalOutOfBagEstimates.count(indexOfSample) == 0) {
                globalOutOfBagEstimates[indexOfSample] = vector<int>(decisionTree->numOutputClasses, 0);
            };
            
            globalOutOfBagEstimates[indexOfSample][predictedOutcomeOfSample] += 1;
        }
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "RandomForest", "updateGlobalOutOfBagEstimates");
        exit(1);
    }  
}
/***********************************************************************/


