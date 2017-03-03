//
//  decisiontree.cpp
//  Mothur
//
//  Created by Sarah Westcott on 10/1/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "decisiontree.hpp"

DecisionTree::DecisionTree(vector< vector<int> >& baseDataSet,
                           vector<int> globalDiscardedFeatureIndices,
                           OptimumFeatureSubsetSelector optimumFeatureSubsetSelector,
                           string treeSplitCriterion,
                           float featureStandardDeviationThreshold)
            : AbstractDecisionTree(baseDataSet,
                                   globalDiscardedFeatureIndices,
                                   optimumFeatureSubsetSelector,
                                   treeSplitCriterion),
            variableImportanceList(numFeatures, 0),
            featureStandardDeviationThreshold(featureStandardDeviationThreshold) {
                
    try {
        m = MothurOut::getInstance();
        createBootStrappedSamples();
        buildDecisionTree();
    }
	catch(exception& e) {
		m->errorOut(e, "DecisionTree", "DecisionTree");
		exit(1);
	} 
}

/***********************************************************************/

int DecisionTree::calcTreeVariableImportanceAndError(int& numCorrect, double& treeErrorRate) {
    try {
        vector< vector<int> > randomlySampledTestData(bootstrappedTestSamples.size(), vector<int>(bootstrappedTestSamples[0].size(), 0));
        
            // TODO: is is possible to further speed up the following O(N^2) by using std::copy?
        for (int i = 0; i < bootstrappedTestSamples.size(); i++) {
            for (int j = 0; j < bootstrappedTestSamples[i].size(); j++) {
                randomlySampledTestData[i][j] = bootstrappedTestSamples[i][j];
            }
        }
        
        for (int i = 0; i < numFeatures; i++) {
            if (m->control_pressed) { return 0; }
            
                // if the index is in globalDiscardedFeatureIndices (i.e, null feature) we don't want to shuffle them
            vector<int>::iterator it = find(globalDiscardedFeatureIndices.begin(), globalDiscardedFeatureIndices.end(), i);
            if (it == globalDiscardedFeatureIndices.end()) {        // NOT FOUND
                // if the standard deviation is very low, we know it's not a good feature at all
                // we can save some time here by discarding that feature
                
                vector<int> featureVector = testSampleFeatureVectors[i];
                if (m->getStandardDeviation(featureVector) > featureStandardDeviationThreshold) {
                    // NOTE: only shuffle the features, never shuffle the output vector
                    // so i = 0 and i will be alwaays <= (numFeatures - 1) as the index at numFeatures will denote
                    // the feature vector
                    randomlyShuffleAttribute(bootstrappedTestSamples, i, i - 1, randomlySampledTestData);

                    int numCorrectAfterShuffle = 0;
                    for (int j = 0; j < randomlySampledTestData.size(); j++) {
                        if (m->control_pressed) {return 0; }
                        
                        vector<int> shuffledSample = randomlySampledTestData[j];
                        int actualSampleOutputClass = shuffledSample[numFeatures];
                        int predictedSampleOutputClass = evaluateSample(shuffledSample);
                        if (actualSampleOutputClass == predictedSampleOutputClass) { numCorrectAfterShuffle++; }
                    }
                    variableImportanceList[i] += (numCorrect - numCorrectAfterShuffle);
                }
            }
        }
        
        // TODO: do we need to save the variableRanks in the DecisionTree, do we need it later?
        vector< pair<int, int> > variableRanks;
        
        for (int i = 0; i < variableImportanceList.size(); i++) {
            if (m->control_pressed) {return 0; }
            if (variableImportanceList[i] > 0) {
                // TODO: is there a way to optimize the follow line's code?
                pair<int, int> variableRank(0, 0);
                variableRank.first = i;
                variableRank.second = variableImportanceList[i];
                variableRanks.push_back(variableRank);
            }
        }
        VariableRankDescendingSorter variableRankDescendingSorter;
        sort(variableRanks.begin(), variableRanks.end(), variableRankDescendingSorter);
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "DecisionTree", "calcTreeVariableImportanceAndError");
		exit(1);
	} 

}
/***********************************************************************/

// TODO: there must be a way to optimize this function
int DecisionTree::evaluateSample(vector<int> testSample) {
    try {
        RFTreeNode *node = rootNode;
        while (true) {
            if (m->control_pressed) { return 0; }
            
            if (node->checkIsLeaf()) { return node->getOutputClass(); }
            
            int sampleSplitFeatureValue = testSample[node->getSplitFeatureIndex()];
            if (sampleSplitFeatureValue < node->getSplitFeatureValue()) { node = node->getLeftChildNode(); }
            else { node = node->getRightChildNode(); } 
        }
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "DecisionTree", "evaluateSample");
		exit(1);
	} 

}
/***********************************************************************/

int DecisionTree::calcTreeErrorRate(int& numCorrect, double& treeErrorRate){
    numCorrect = 0;
    try {
        for (int i = 0; i < bootstrappedTestSamples.size(); i++) {
             if (m->control_pressed) {return 0; }
            
            vector<int> testSample = bootstrappedTestSamples[i];
            int testSampleIndex = bootstrappedTestSampleIndices[i];
            
            int actualSampleOutputClass = testSample[numFeatures];
            int predictedSampleOutputClass = evaluateSample(testSample);
            
            if (actualSampleOutputClass == predictedSampleOutputClass) { numCorrect++; } 
            
            outOfBagEstimates[testSampleIndex] = predictedSampleOutputClass;
        }
        
        treeErrorRate = 1 - ((double)numCorrect / (double)bootstrappedTestSamples.size());   
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "DecisionTree", "calcTreeErrorRate");
		exit(1);
	} 
}

/***********************************************************************/
// TODO: optimize the algo, instead of transposing two time, we can extarct the feature,
// shuffle it and then re-insert in the original place, thus iproving runnting time
//This function randomize abundances for a given OTU/feature.

void DecisionTree::randomlyShuffleAttribute(const vector< vector<int> >& samples,
                               const int featureIndex,
                               const int prevFeatureIndex,
                               vector< vector<int> >& shuffledSample) {
    try {
        // NOTE: we need (numFeatures + 1) featureVecotors, the last extra vector is actually outputVector
        
        // restore previously shuffled feature
        if (prevFeatureIndex > -1) {
            for (int j = 0; j < samples.size(); j++) {
                if (m->control_pressed) { return; }
                shuffledSample[j][prevFeatureIndex] = samples[j][prevFeatureIndex];
            }
        }
        
        // now do the shuffling
        vector<int> featureVectors(samples.size(), 0);
        for (int j = 0; j < samples.size(); j++) {
            if (m->control_pressed) { return; }
            featureVectors[j] = samples[j][featureIndex];
        }
        m->mothurRandomShuffle(featureVectors);
        for (int j = 0; j < samples.size(); j++) {
            if (m->control_pressed) { return; }
            shuffledSample[j][featureIndex] = featureVectors[j];
        }
    }
	catch(exception& e) {
        m->errorOut(e, "DecisionTree", "randomlyShuffleAttribute");
		exit(1);
	}
    
}

/***********************************************************************/

int DecisionTree::purgeTreeNodesDataRecursively(RFTreeNode* treeNode) {
    try {
        treeNode->bootstrappedTrainingSamples.clear();
        treeNode->bootstrappedFeatureVectors.clear();
        treeNode->bootstrappedOutputVector.clear();
        treeNode->localDiscardedFeatureIndices.clear();
        treeNode->globalDiscardedFeatureIndices.clear();
        
        if (treeNode->leftChildNode != NULL) { purgeTreeNodesDataRecursively(treeNode->leftChildNode); }
        if (treeNode->rightChildNode != NULL) { purgeTreeNodesDataRecursively(treeNode->rightChildNode); }
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "DecisionTree", "purgeTreeNodesDataRecursively");
		exit(1);
	} 
}
/***********************************************************************/

void DecisionTree::buildDecisionTree(){
    try {
    
        int generation = 0;
        rootNode = new RFTreeNode(bootstrappedTrainingSamples, globalDiscardedFeatureIndices, numFeatures, numSamples, numOutputClasses, generation, nodeIdCount, featureStandardDeviationThreshold);
        nodeIdCount++;
        
        splitRecursively(rootNode);
        
        }
	catch(exception& e) {
		m->errorOut(e, "DecisionTree", "buildDecisionTree");
		exit(1);
	} 
}

/***********************************************************************/

int DecisionTree::splitRecursively(RFTreeNode* rootNode) {
    try {
       
        if (rootNode->getNumSamples() < 2){
            rootNode->setIsLeaf(true);
            rootNode->setOutputClass(rootNode->getBootstrappedTrainingSamples()[0][rootNode->getNumFeatures()]);
            return 0;
        }
        
        int classifiedOutputClass;
        bool isAlreadyClassified = checkIfAlreadyClassified(rootNode, classifiedOutputClass);    
        if (isAlreadyClassified == true){
            rootNode->setIsLeaf(true);
            rootNode->setOutputClass(classifiedOutputClass);
            return 0;
        }
        if (m->control_pressed) { return 0; }
        vector<int> featureSubsetIndices = selectFeatureSubsetRandomly(globalDiscardedFeatureIndices, rootNode->getLocalDiscardedFeatureIndices());
        
            // TODO: need to check if the value is actually copied correctly
        rootNode->setFeatureSubsetIndices(featureSubsetIndices);
        if (m->control_pressed) { return 0; }
      
        findAndUpdateBestFeatureToSplitOn(rootNode);
        
        // update rootNode outputClass, this is needed for pruning
        // this is only for internal nodes
        updateOutputClassOfNode(rootNode);
        
        if (m->control_pressed) { return 0; }
        
        vector< vector<int> > leftChildSamples;
        vector< vector<int> > rightChildSamples;
        getSplitPopulation(rootNode, leftChildSamples, rightChildSamples);
        
        if (m->control_pressed) { return 0; }
        
        // TODO: need to write code to clear this memory
        RFTreeNode* leftChildNode = new RFTreeNode(leftChildSamples, globalDiscardedFeatureIndices, numFeatures, (int)leftChildSamples.size(), numOutputClasses, rootNode->getGeneration() + 1, nodeIdCount, featureStandardDeviationThreshold);
        nodeIdCount++;
        RFTreeNode* rightChildNode = new RFTreeNode(rightChildSamples, globalDiscardedFeatureIndices, numFeatures, (int)rightChildSamples.size(), numOutputClasses, rootNode->getGeneration() + 1, nodeIdCount, featureStandardDeviationThreshold);
        nodeIdCount++;
        
        rootNode->setLeftChildNode(leftChildNode);
        leftChildNode->setParentNode(rootNode);
        
        rootNode->setRightChildNode(rightChildNode);
        rightChildNode->setParentNode(rootNode);
        
        // TODO: This recursive split can be parrallelized later
        splitRecursively(leftChildNode);
        if (m->control_pressed) { return 0; }
        
        splitRecursively(rightChildNode);
        return 0;
        
    }
	catch(exception& e) {
		m->errorOut(e, "DecisionTree", "splitRecursively");
		exit(1);
	} 
}
/***********************************************************************/

int DecisionTree::findAndUpdateBestFeatureToSplitOn(RFTreeNode* node){
    try {

        vector< vector<int> > bootstrappedFeatureVectors = node->getBootstrappedFeatureVectors();
        if (m->control_pressed) { return 0; }
        vector<int> bootstrappedOutputVector = node->getBootstrappedOutputVector();
        if (m->control_pressed) { return 0; }
        vector<int> featureSubsetIndices = node->getFeatureSubsetIndices();
        if (m->control_pressed) { return 0; }
        
        vector<double> featureSubsetEntropies;
        vector<int> featureSubsetSplitValues;
        vector<double> featureSubsetIntrinsicValues;
        vector<double> featureSubsetGainRatios;
        
        for (int i = 0; i < featureSubsetIndices.size(); i++) {
            if (m->control_pressed) { return 0; }
            
            int tryIndex = featureSubsetIndices[i];
                       
            double featureMinEntropy;
            int featureSplitValue;
            double featureIntrinsicValue;
            
            getMinEntropyOfFeature(bootstrappedFeatureVectors[tryIndex], bootstrappedOutputVector, featureMinEntropy, featureSplitValue, featureIntrinsicValue);
            if (m->control_pressed) { return 0; }
            
            featureSubsetEntropies.push_back(featureMinEntropy);
            featureSubsetSplitValues.push_back(featureSplitValue);
            featureSubsetIntrinsicValues.push_back(featureIntrinsicValue);
            
            double featureInformationGain = node->getOwnEntropy() - featureMinEntropy;
            double featureGainRatio = (double)featureInformationGain / (double)featureIntrinsicValue;
            featureSubsetGainRatios.push_back(featureGainRatio);
            
        }
        
        vector<double>::iterator minEntropyIterator = min_element(featureSubsetEntropies.begin(), featureSubsetEntropies.end());
        vector<double>::iterator maxGainRatioIterator = max_element(featureSubsetGainRatios.begin(), featureSubsetGainRatios.end());
        double featureMinEntropy = *minEntropyIterator;
        
        // TODO: kept the following line as future reference, can be useful
        // double featureMaxGainRatio = *maxGainRatioIterator;
        
        double bestFeatureSplitEntropy = featureMinEntropy;
        int bestFeatureToSplitOnIndex = -1;
        if (treeSplitCriterion == "gainratio"){
            bestFeatureToSplitOnIndex = (int)(maxGainRatioIterator - featureSubsetGainRatios.begin());
            // if using 'gainRatio' measure, then featureMinEntropy must be re-updated, as the index
            // for 'featureMaxGainRatio' would be different
            bestFeatureSplitEntropy = featureSubsetEntropies[bestFeatureToSplitOnIndex];
        } else  if ( treeSplitCriterion == "infogain"){
            bestFeatureToSplitOnIndex = (int)(minEntropyIterator - featureSubsetEntropies.begin());
        } else {
                // TODO: we need an abort mechanism here
        }
        
            // TODO: is the following line needed? kept is as future reference
        // splitInformationGain = node.ownEntropy - node.splitFeatureEntropy
        
        int bestFeatureSplitValue = featureSubsetSplitValues[bestFeatureToSplitOnIndex];
        
        node->setSplitFeatureIndex(featureSubsetIndices[bestFeatureToSplitOnIndex]);
        node->setSplitFeatureValue(bestFeatureSplitValue);
        node->setSplitFeatureEntropy(bestFeatureSplitEntropy);
            // TODO: kept the following line as future reference
        // node.splitInformationGain = splitInformationGain
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "DecisionTree", "findAndUpdateBestFeatureToSplitOn");
		exit(1);
	} 
}
/***********************************************************************/
vector<int> DecisionTree::selectFeatureSubsetRandomly(vector<int> globalDiscardedFeatureIndices, vector<int> localDiscardedFeatureIndices){
    try {

        vector<int> featureSubsetIndices;
        
        vector<int> combinedDiscardedFeatureIndices;
        combinedDiscardedFeatureIndices.insert(combinedDiscardedFeatureIndices.end(), globalDiscardedFeatureIndices.begin(), globalDiscardedFeatureIndices.end());
        combinedDiscardedFeatureIndices.insert(combinedDiscardedFeatureIndices.end(), localDiscardedFeatureIndices.begin(), localDiscardedFeatureIndices.end());
        
        sort(combinedDiscardedFeatureIndices.begin(), combinedDiscardedFeatureIndices.end());
        
        int numberOfRemainingSuitableFeatures = (int)(numFeatures - combinedDiscardedFeatureIndices.size());
        int currentFeatureSubsetSize = numberOfRemainingSuitableFeatures < optimumFeatureSubsetSize ? numberOfRemainingSuitableFeatures : optimumFeatureSubsetSize;
        
        while (featureSubsetIndices.size() < currentFeatureSubsetSize) {
            
            if (m->control_pressed) { return featureSubsetIndices; }
            
            int randomIndex = m->getRandomIndex(numFeatures-1);
            vector<int>::iterator it = find(featureSubsetIndices.begin(), featureSubsetIndices.end(), randomIndex);
            if (it == featureSubsetIndices.end()){    // NOT FOUND
                vector<int>::iterator it2 = find(combinedDiscardedFeatureIndices.begin(), combinedDiscardedFeatureIndices.end(), randomIndex);
                if (it2 == combinedDiscardedFeatureIndices.end()){  // NOT FOUND AGAIN
                    featureSubsetIndices.push_back(randomIndex);
                }
            }
        }
        sort(featureSubsetIndices.begin(), featureSubsetIndices.end());
        
        //#ifdef DEBUG_LEVEL_3
        //    PRINT_VAR(featureSubsetIndices);
        //#endif
        
        return featureSubsetIndices;
    }
	catch(exception& e) {
		m->errorOut(e, "DecisionTree", "selectFeatureSubsetRandomly");
		exit(1);
	} 
}
/***********************************************************************/

// TODO: printTree() needs a check if correct
int DecisionTree::printTree(RFTreeNode* treeNode, string caption){
    try { 
        string tabs = "";
        for (int i = 0; i < treeNode->getGeneration(); i++) { tabs += "|--"; }
        //    for (int i = 0; i < treeNode->getGeneration() - 1; i++) { tabs += "|  "; }
        //    if (treeNode->getGeneration() != 0) { tabs += "|--"; }
        
        if (treeNode != NULL && treeNode->checkIsLeaf() == false){
            m->mothurOut(tabs + caption + " [ gen: " + toString(treeNode->getGeneration()) + " , id: " + toString(treeNode->nodeId) + " ] ( " + toString(treeNode->getSplitFeatureValue()) + " < X" + toString(treeNode->getSplitFeatureIndex()) + " ) ( predicted: " + toString(treeNode->outputClass) + " , misclassified: " + toString(treeNode->testSampleMisclassificationCount) + " )\n");
            
            printTree(treeNode->getLeftChildNode(), "left ");
            printTree(treeNode->getRightChildNode(), "right");
        }else {
            m->mothurOut(tabs + caption + " [ gen: " + toString(treeNode->getGeneration()) + " , id: " + toString(treeNode->nodeId) + " ] ( classified to: " + toString(treeNode->getOutputClass()) + ", samples: " + toString(treeNode->getNumSamples()) + " , misclassified: " + toString(treeNode->testSampleMisclassificationCount) + " )\n");
        }
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "DecisionTree", "printTree");
		exit(1);
	} 
}
/***********************************************************************/
void DecisionTree::deleteTreeNodesRecursively(RFTreeNode* treeNode) {
    try {
        if (treeNode == NULL) { return; }
        deleteTreeNodesRecursively(treeNode->leftChildNode);
        deleteTreeNodesRecursively(treeNode->rightChildNode);
        delete treeNode; treeNode = NULL;
    }
	catch(exception& e) {
		m->errorOut(e, "DecisionTree", "deleteTreeNodesRecursively");
		exit(1);
	} 
}
/***********************************************************************/

void DecisionTree::pruneTree(double pruneAggressiveness = 0.9) {
    
    // find out the number of misclassification by each of the nodes
    for (int i = 0; i < bootstrappedTestSamples.size(); i++) {
        if (m->control_pressed) { return; }
        
        vector<int> testSample = bootstrappedTestSamples[i];
        updateMisclassificationCountRecursively(rootNode, testSample);
    }
    
    // do the actual pruning
    pruneRecursively(rootNode, pruneAggressiveness);
}
/***********************************************************************/

void DecisionTree::pruneRecursively(RFTreeNode* treeNode, double pruneAggressiveness){
    
    if (treeNode != NULL && treeNode->checkIsLeaf() == false) {
        if (m->control_pressed) { return; }
        
        pruneRecursively(treeNode->leftChildNode, pruneAggressiveness);
        pruneRecursively(treeNode->rightChildNode, pruneAggressiveness);
        
        int subTreeMisclassificationCount = treeNode->leftChildNode->getTestSampleMisclassificationCount() + treeNode->rightChildNode->getTestSampleMisclassificationCount();
        int ownMisclassificationCount = treeNode->getTestSampleMisclassificationCount();
        
        if (subTreeMisclassificationCount * pruneAggressiveness > ownMisclassificationCount) {
                // TODO: need to check the effect of these two delete calls
            delete treeNode->leftChildNode;
            treeNode->leftChildNode = NULL;
            
            delete treeNode->rightChildNode;
            treeNode->rightChildNode = NULL;
            
            treeNode->isLeaf = true;
        }
        
    }
}
/***********************************************************************/

void DecisionTree::updateMisclassificationCountRecursively(RFTreeNode* treeNode, vector<int> testSample) {
    
    int actualSampleOutputClass = testSample[numFeatures];
    int nodePredictedOutputClass = treeNode->outputClass;
    
    if (actualSampleOutputClass != nodePredictedOutputClass) {
        treeNode->testSampleMisclassificationCount++;
        map<int, int>::iterator it = nodeMisclassificationCounts.find(treeNode->nodeId);
        if (it == nodeMisclassificationCounts.end()) {  // NOT FOUND
            nodeMisclassificationCounts[treeNode->nodeId] = 0;
        }
        nodeMisclassificationCounts[treeNode->nodeId]++;
    }
    
    if (treeNode->checkIsLeaf() == false) { // NOT A LEAF
        int sampleSplitFeatureValue = testSample[treeNode->splitFeatureIndex];
        if (sampleSplitFeatureValue < treeNode->splitFeatureValue) {
            updateMisclassificationCountRecursively(treeNode->leftChildNode, testSample);
        } else {
            updateMisclassificationCountRecursively(treeNode->rightChildNode, testSample);
        }
    }
}

/***********************************************************************/

void DecisionTree::updateOutputClassOfNode(RFTreeNode* treeNode) {
    vector<int> counts(numOutputClasses, 0);
    for (int i = 0; i < treeNode->bootstrappedOutputVector.size(); i++) {
        int bootstrappedOutput = treeNode->bootstrappedOutputVector[i];
        counts[bootstrappedOutput]++;
    }

    vector<int>::iterator majorityVotedOutputClassCountIterator = max_element(counts.begin(), counts.end());
    int majorityVotedOutputClassCount = *majorityVotedOutputClassCountIterator;
    vector<int>::iterator it = find(counts.begin(), counts.end(), majorityVotedOutputClassCount);
    int majorityVotedOutputClass = (int)(it - counts.begin());
    treeNode->setOutputClass(majorityVotedOutputClass);

}
/***********************************************************************/


