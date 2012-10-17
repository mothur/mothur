//
//  decisiontree.cpp
//  Mothur
//
//  Created by Sarah Westcott on 10/1/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "decisiontree.hpp"

DecisionTree::DecisionTree(vector< vector<int> > baseDataSet,
             vector<int> globalDiscardedFeatureIndices,
             OptimumFeatureSubsetSelector optimumFeatureSubsetSelector,
             string treeSplitCriterion) : AbstractDecisionTree(baseDataSet,
                       globalDiscardedFeatureIndices,
                       optimumFeatureSubsetSelector,
                       treeSplitCriterion), variableImportanceList(numFeatures, 0){
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

int DecisionTree::calcTreeVariableImportanceAndError() {
    try {
        
        int numCorrect;
        double treeErrorRate;
        calcTreeErrorRate(numCorrect, treeErrorRate);
        
        if (m->control_pressed) {return 0; }
                
        for (int i = 0; i < numFeatures; i++) {
            if (m->control_pressed) {return 0; }
            // NOTE: only shuffle the features, never shuffle the output vector
            // so i = 0 and i will be alwaays <= (numFeatures - 1) as the index at numFeatures will denote
            // the feature vector
            vector< vector<int> > randomlySampledTestData = randomlyShuffleAttribute(bootstrappedTestSamples, i);
            
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
        
        // TODO: do we need to save the variableRanks in the DecisionTree, do we need it later?
        vector< vector<int> > variableRanks;
        for (int i = 0; i < variableImportanceList.size(); i++) {
            if (m->control_pressed) {return 0; }
            if (variableImportanceList[i] > 0) {
                // TODO: is there a way to optimize the follow line's code?
                vector<int> variableRank(2, 0);
                variableRank[0] = i; variableRank[1] = variableImportanceList[i];
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
            if (m->control_pressed) {return 0; }
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
    try {
        numCorrect = 0;
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
vector< vector<int> > DecisionTree::randomlyShuffleAttribute(vector< vector<int> > samples, int featureIndex) {
    try {
        // NOTE: we need (numFeatures + 1) featureVecotors, the last extra vector is actually outputVector
        vector< vector<int> > shuffledSample = samples;
        vector<int> featureVectors(samples.size(), 0);
        
        for (int j = 0; j < samples.size(); j++) {
            if (m->control_pressed) { return shuffledSample; }
            featureVectors[j] = samples[j][featureIndex];
        }
        
        random_shuffle(featureVectors.begin(), featureVectors.end());

        for (int j = 0; j < samples.size(); j++) {
            if (m->control_pressed) {return shuffledSample; }
            shuffledSample[j][featureIndex] = featureVectors[j];
        }
        
        return shuffledSample;
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
    rootNode = new RFTreeNode(bootstrappedTrainingSamples, globalDiscardedFeatureIndices, numFeatures, numSamples, numOutputClasses, generation);
    
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
        if (m->control_pressed) {return 0;}
        vector<int> featureSubsetIndices = selectFeatureSubsetRandomly(globalDiscardedFeatureIndices, rootNode->getLocalDiscardedFeatureIndices());
        rootNode->setFeatureSubsetIndices(featureSubsetIndices);
        if (m->control_pressed) {return 0;}
      
        findAndUpdateBestFeatureToSplitOn(rootNode);
        
        if (m->control_pressed) {return 0;}
        
        vector< vector<int> > leftChildSamples;
        vector< vector<int> > rightChildSamples;
        getSplitPopulation(rootNode, leftChildSamples, rightChildSamples);
        
        if (m->control_pressed) {return 0;}
        
        // TODO: need to write code to clear this memory
        RFTreeNode* leftChildNode = new RFTreeNode(leftChildSamples, globalDiscardedFeatureIndices, numFeatures, (int)leftChildSamples.size(), numOutputClasses, rootNode->getGeneration() + 1);
        RFTreeNode* rightChildNode = new RFTreeNode(rightChildSamples, globalDiscardedFeatureIndices, numFeatures, (int)rightChildSamples.size(), numOutputClasses, rootNode->getGeneration() + 1);
        
        rootNode->setLeftChildNode(leftChildNode);
        leftChildNode->setParentNode(rootNode);
        
        rootNode->setRightChildNode(rightChildNode);
        rightChildNode->setParentNode(rootNode);
        
        // TODO: This recursive split can be parrallelized later
        splitRecursively(leftChildNode);
        if (m->control_pressed) {return 0;}
        
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
        if (m->control_pressed) {return 0;}
        vector<int> bootstrappedOutputVector = node->getBootstrappedOutputVector();
        if (m->control_pressed) {return 0;}
        vector<int> featureSubsetIndices = node->getFeatureSubsetIndices();
        if (m->control_pressed) {return 0;}
        
        vector<double> featureSubsetEntropies;
        vector<int> featureSubsetSplitValues;
        vector<double> featureSubsetIntrinsicValues;
        vector<double> featureSubsetGainRatios;
        
        for (int i = 0; i < featureSubsetIndices.size(); i++) {
            if (m->control_pressed) {return 0;}
            
            int tryIndex = featureSubsetIndices[i];
                       
            double featureMinEntropy;
            int featureSplitValue;
            double featureIntrinsicValue;
            
            getMinEntropyOfFeature(bootstrappedFeatureVectors[tryIndex], bootstrappedOutputVector, featureMinEntropy, featureSplitValue, featureIntrinsicValue);
            if (m->control_pressed) {return 0;}
            
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
        //double featureMaxGainRatio = *maxGainRatioIterator;
        
        double bestFeatureSplitEntropy = featureMinEntropy;
        int bestFeatureToSplitOnIndex = -1;
        if (treeSplitCriterion == "gainRatio"){ 
            bestFeatureToSplitOnIndex = (int)(maxGainRatioIterator - featureSubsetGainRatios.begin());
            // if using 'gainRatio' measure, then featureMinEntropy must be re-updated, as the index
            // for 'featureMaxGainRatio' would be different
            bestFeatureSplitEntropy = featureSubsetEntropies[bestFeatureToSplitOnIndex];
        }
        else { bestFeatureToSplitOnIndex = (int)(minEntropyIterator - featureSubsetEntropies.begin()); }
        
        int bestFeatureSplitValue = featureSubsetSplitValues[bestFeatureToSplitOnIndex];
        
        node->setSplitFeatureIndex(featureSubsetIndices[bestFeatureToSplitOnIndex]);
        node->setSplitFeatureValue(bestFeatureSplitValue);
        node->setSplitFeatureEntropy(bestFeatureSplitEntropy);
        
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
            
            // TODO: optimize rand() call here
            int randomIndex = rand() % numFeatures;
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
        for (int i = 0; i < treeNode->getGeneration(); i++) { tabs += "   "; }
        //    for (int i = 0; i < treeNode->getGeneration() - 1; i++) { tabs += "|  "; }
        //    if (treeNode->getGeneration() != 0) { tabs += "|--"; }
        
        if (treeNode != NULL && treeNode->checkIsLeaf() == false){
            m->mothurOut(tabs + caption + " [ gen: " + toString(treeNode->getGeneration()) + " ] ( " + toString(treeNode->getSplitFeatureValue()) + " < X" + toString(treeNode->getSplitFeatureIndex()) +" )\n");
            
            printTree(treeNode->getLeftChildNode(), "leftChild");
            printTree(treeNode->getRightChildNode(), "rightChild");
        }else {
            m->mothurOut(tabs + caption + " [ gen: " + toString(treeNode->getGeneration()) + " ] ( classified to: " + toString(treeNode->getOutputClass()) + ", samples: " + toString(treeNode->getNumSamples()) + " )\n");
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
        delete treeNode;
    }
	catch(exception& e) {
		m->errorOut(e, "DecisionTree", "deleteTreeNodesRecursively");
		exit(1);
	} 
}
/***********************************************************************/

