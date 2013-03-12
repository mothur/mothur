//
//  randomforest.cpp
//  Mothur
//
//  Created by Sarah Westcott on 10/2/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "randomforest.hpp" 

/***********************************************************************/

RandomForest::RandomForest(const vector <vector<int> > dataSet,const int numDecisionTrees,
             const string treeSplitCriterion = "informationGain") : Forest(dataSet, numDecisionTrees, treeSplitCriterion) {
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
    
    // TODO: need to add try/catch operators to fix this
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
// DONE
int RandomForest::populateDecisionTrees() {
    try {
        
        for (int i = 0; i < numDecisionTrees; i++) {
            if (m->control_pressed) { return 0; }
            if (((i+1) % 10) == 0) {  m->mothurOut("Creating " + toString(i+1) + " (th) Decision tree\n");  }
            // TODO: need to first fix if we are going to use pointer based system or anything else
            DecisionTree* decisionTree = new DecisionTree(dataSet, globalDiscardedFeatureIndices, OptimumFeatureSubsetSelector("log2"), treeSplitCriterion);
            decisionTree->calcTreeVariableImportanceAndError();
            if (m->control_pressed) { return 0; }
            updateGlobalOutOfBagEstimates(decisionTree);
            if (m->control_pressed) { return 0; }
            decisionTree->purgeDataSetsFromTree();
            if (m->control_pressed) { return 0; }
            decisionTrees.push_back(decisionTree);
        }
        
        if (m->debug) {
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


