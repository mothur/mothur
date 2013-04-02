//
//  abstractdecisiontree.cpp
//  Mothur
//
//  Created by Sarah Westcott on 10/1/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "abstractdecisiontree.hpp"

/**************************************************************************************************/

AbstractDecisionTree::AbstractDecisionTree(vector<vector<int> > baseDataSet, 
                                         vector<int> globalDiscardedFeatureIndices,
                                         OptimumFeatureSubsetSelector optimumFeatureSubsetSelector, 
                                         string treeSplitCriterion)

                    : baseDataSet(baseDataSet),
                    numSamples((int)baseDataSet.size()),
                    numFeatures((int)(baseDataSet[0].size() - 1)),
                    numOutputClasses(0),
                    rootNode(NULL),
                    nodeIdCount(0),
                    globalDiscardedFeatureIndices(globalDiscardedFeatureIndices),
                    optimumFeatureSubsetSize(optimumFeatureSubsetSelector.getOptimumFeatureSubsetSize(numFeatures)),
                    treeSplitCriterion(treeSplitCriterion) {

    try {
        // TODO: istead of calculating this for every DecisionTree
        // clacualte this once in the RandomForest class and pass the values
        m = MothurOut::getInstance();
        for (int i = 0;  i < numSamples; i++) {
            if (m->control_pressed) { break; }
            int outcome = baseDataSet[i][numFeatures];
            vector<int>::iterator it = find(outputClasses.begin(), outputClasses.end(), outcome);
            if (it == outputClasses.end()){       // find() will return classes.end() if the element is not found
                outputClasses.push_back(outcome);
                numOutputClasses++;
            }
        }
        
        if (m->debug) {
            //m->mothurOut("outputClasses = " + toStringVectorInt(outputClasses));
            m->mothurOut("numOutputClasses = " + toString(numOutputClasses) + '\n');
        }

    }
	catch(exception& e) {
		m->errorOut(e, "AbstractDecisionTree", "AbstractDecisionTree");
		exit(1);
	} 
}
/**************************************************************************************************/
int AbstractDecisionTree::createBootStrappedSamples(){
    try {    
        vector<bool> isInTrainingSamples(numSamples, false);
        
        for (int i = 0; i < numSamples; i++) {
            if (m->control_pressed) { return 0; }
            // TODO: optimize the rand() function call + double check if it's working properly
            int randomIndex = rand() % numSamples;
            bootstrappedTrainingSamples.push_back(baseDataSet[randomIndex]);
            isInTrainingSamples[randomIndex] = true;
        }
        
        for (int i = 0; i < numSamples; i++) {
            if (m->control_pressed) { return 0; }
            if (isInTrainingSamples[i]){ bootstrappedTrainingSampleIndices.push_back(i); }
            else{
                bootstrappedTestSamples.push_back(baseDataSet[i]);
                bootstrappedTestSampleIndices.push_back(i);
            }
        }
        
            // do the transpose of Test Samples
        for (int i = 0; i < bootstrappedTestSamples[0].size(); i++) {
            if (m->control_pressed) { return 0; }
            
            vector<int> tmpFeatureVector(bootstrappedTestSamples.size(), 0);
            for (int j = 0; j < bootstrappedTestSamples.size(); j++) {
                if (m->control_pressed) { return 0; }
                
                tmpFeatureVector[j] = bootstrappedTestSamples[j][i];
            }
            testSampleFeatureVectors.push_back(tmpFeatureVector);
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "AbstractDecisionTree", "createBootStrappedSamples");
		exit(1);
	} 
}
/**************************************************************************************************/
int AbstractDecisionTree::getMinEntropyOfFeature(vector<int> featureVector,
                                                 vector<int> outputVector,
                                                 double& minEntropy,
                                                 int& featureSplitValue,
                                                 double& intrinsicValue){
    try {

        vector< pair<int, int> > featureOutputPair(featureVector.size(), pair<int, int>(0, 0));
        
        for (int i = 0; i < featureVector.size(); i++) { 
            if (m->control_pressed) { return 0; }
            
            featureOutputPair[i].first = featureVector[i];
            featureOutputPair[i].second = outputVector[i];
        }
        // TODO: using default behavior to sort(), need to specify the comparator for added safety and compiler portability,
        
        IntPairVectorSorter intPairVectorSorter;
        sort(featureOutputPair.begin(), featureOutputPair.end(), intPairVectorSorter);
        
        vector<int> splitPoints;
        vector<int> uniqueFeatureValues(1, featureOutputPair[0].first);
        
        for (int i = 0; i < featureOutputPair.size(); i++) {

            if (m->control_pressed) { return 0; }
            int featureValue = featureOutputPair[i].first;

            vector<int>::iterator it = find(uniqueFeatureValues.begin(), uniqueFeatureValues.end(), featureValue);
            if (it == uniqueFeatureValues.end()){                 // NOT FOUND
                uniqueFeatureValues.push_back(featureValue);
                splitPoints.push_back(i);
            }
        }
        

        
        int bestSplitIndex = -1;
        if (splitPoints.size() == 0){
            // TODO: trying out C++'s infitinity, don't know if this will work properly
            // TODO: check the caller function of this function, there check the value if minEntropy and comapre to inf
            // so that no wrong calculation is done
            minEntropy = numeric_limits<double>::infinity();                          // OUTPUT
            intrinsicValue = numeric_limits<double>::infinity();                      // OUTPUT
            featureSplitValue = -1;                                                   // OUTPUT
        }else{
            getBestSplitAndMinEntropy(featureOutputPair, splitPoints, minEntropy, bestSplitIndex, intrinsicValue);  // OUTPUT
            featureSplitValue = featureOutputPair[splitPoints[bestSplitIndex]].first;    // OUTPUT
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "AbstractDecisionTree", "getMinEntropyOfFeature");
		exit(1);
	} 
}
/**************************************************************************************************/
double AbstractDecisionTree::calcIntrinsicValue(int numLessThanValueAtSplitPoint, int numGreaterThanValueAtSplitPoint, int numSamples) {
    try {
        double upperSplitEntropy = 0.0, lowerSplitEntropy = 0.0;
        if (numLessThanValueAtSplitPoint > 0) {
            upperSplitEntropy = numLessThanValueAtSplitPoint * log2((double) numLessThanValueAtSplitPoint / (double) numSamples);
        }
        
        if (numGreaterThanValueAtSplitPoint > 0) {
            lowerSplitEntropy = numGreaterThanValueAtSplitPoint * log2((double) numGreaterThanValueAtSplitPoint / (double) numSamples);
        }
        
        double intrinsicValue = - ((double)(upperSplitEntropy + lowerSplitEntropy) / (double)numSamples);
        return intrinsicValue;
    }
	catch(exception& e) {
		m->errorOut(e, "AbstractDecisionTree", "calcIntrinsicValue");
		exit(1);
	} 
}
/**************************************************************************************************/

int AbstractDecisionTree::getBestSplitAndMinEntropy(vector< pair<int, int> > featureOutputPairs, vector<int> splitPoints,
                                                    double& minEntropy, int& minEntropyIndex, double& relatedIntrinsicValue){
    try {
        
        int numSamples = (int)featureOutputPairs.size();
        vector<double> entropies;
        vector<double> intrinsicValues;
        
        for (int i = 0; i < splitPoints.size(); i++) {
            if (m->control_pressed) { return 0; }
            int index = splitPoints[i];
            int valueAtSplitPoint = featureOutputPairs[index].first;

            int numLessThanValueAtSplitPoint = 0;
            int numGreaterThanValueAtSplitPoint = 0;
            
            for (int j = 0; j < featureOutputPairs.size(); j++) {
                if (m->control_pressed) { return 0; }
                pair<int, int> record = featureOutputPairs[j];
                if (record.first < valueAtSplitPoint){ numLessThanValueAtSplitPoint++; }
                else{ numGreaterThanValueAtSplitPoint++; }
            }
            
            double upperEntropyOfSplit = calcSplitEntropy(featureOutputPairs, index, numOutputClasses, true);
            double lowerEntropyOfSplit = calcSplitEntropy(featureOutputPairs, index, numOutputClasses, false);
            
            double totalEntropy = (numLessThanValueAtSplitPoint * upperEntropyOfSplit + numGreaterThanValueAtSplitPoint * lowerEntropyOfSplit) / (double)numSamples;
            double intrinsicValue = calcIntrinsicValue(numLessThanValueAtSplitPoint, numGreaterThanValueAtSplitPoint, numSamples);
            entropies.push_back(totalEntropy);
            intrinsicValues.push_back(intrinsicValue);
            
        }
                
        // set output values
        vector<double>::iterator it = min_element(entropies.begin(), entropies.end());
        minEntropy = *it;                                                         // OUTPUT
        minEntropyIndex = (int)(it - entropies.begin());                          // OUTPUT
        relatedIntrinsicValue = intrinsicValues[minEntropyIndex];                 // OUTPUT
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "AbstractDecisionTree", "getBestSplitAndMinEntropy");
		exit(1);
	} 
}
/**************************************************************************************************/

double AbstractDecisionTree::calcSplitEntropy(vector< pair<int, int> > featureOutputPairs, int splitIndex, int numOutputClasses, bool isUpperSplit = true) {
    try {
        vector<int> classCounts(numOutputClasses, 0);
        
        if (isUpperSplit) { 
            for (int i = 0; i < splitIndex; i++) {
                if (m->control_pressed) { return 0; }
                classCounts[featureOutputPairs[i].second]++;
            }
        } else {
            for (int i = splitIndex; i < featureOutputPairs.size(); i++) { 
                if (m->control_pressed) { return 0; }
                classCounts[featureOutputPairs[i].second]++;
            }
        }
        
        int totalClassCounts = accumulate(classCounts.begin(), classCounts.end(), 0);
        
        double splitEntropy = 0.0;
        
        for (int i = 0; i < classCounts.size(); i++) {
            if (m->control_pressed) { return 0; }
            if (classCounts[i] == 0) { continue; }
            double probability = (double) classCounts[i] / (double) totalClassCounts;
            splitEntropy += -(probability * log2(probability));
        }
        
        return splitEntropy;
    }
	catch(exception& e) {
		m->errorOut(e, "AbstractDecisionTree", "calcSplitEntropy");
		exit(1);
	} 
}

/**************************************************************************************************/

int AbstractDecisionTree::getSplitPopulation(RFTreeNode* node, vector< vector<int> >& leftChildSamples, vector< vector<int> >& rightChildSamples){    
    try {
        // TODO: there is a possibility of optimization if we can recycle the samples in each nodes
        // we just need to pointers to the samples i.e. vector<int> and use it everywhere and not create the sample 
        // sample over and over again
        // we need to make this const so that it is not modified by all the function calling
        // currently purgeTreeNodesDataRecursively() is used for the same purpose, but this can be avoided altogher
        // if re-using the same data over the classes
        
        int splitFeatureGlobalIndex = node->getSplitFeatureIndex();
        
        for (int i = 0; i < node->getBootstrappedTrainingSamples().size(); i++) {
            if (m->control_pressed) { return 0; }
            vector<int> sample =  node->getBootstrappedTrainingSamples()[i];
            if (m->control_pressed) { return 0; }
            
            if (sample[splitFeatureGlobalIndex] < node->getSplitFeatureValue()) { leftChildSamples.push_back(sample); }
            else { rightChildSamples.push_back(sample); }
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "AbstractDecisionTree", "getSplitPopulation");
		exit(1);
	} 
}
/**************************************************************************************************/
// TODO: checkIfAlreadyClassified() verify code
// TODO: use bootstrappedOutputVector for easier calculation instead of using getBootstrappedTrainingSamples()
bool AbstractDecisionTree::checkIfAlreadyClassified(RFTreeNode* treeNode, int& outputClass) {
    try {

        vector<int> tempOutputClasses;
        for (int i = 0; i < treeNode->getBootstrappedTrainingSamples().size(); i++) {
            if (m->control_pressed) { return 0; }
            int sampleOutputClass = treeNode->getBootstrappedTrainingSamples()[i][numFeatures];
            vector<int>::iterator it = find(tempOutputClasses.begin(), tempOutputClasses.end(), sampleOutputClass);
            if (it == tempOutputClasses.end()) {               // NOT FOUND
                tempOutputClasses.push_back(sampleOutputClass);
            }
        }
        
        if (tempOutputClasses.size() < 2) { outputClass = tempOutputClasses[0]; return true; }
        else { outputClass = -1; return false; }
        
    }
	catch(exception& e) {
		m->errorOut(e, "AbstractDecisionTree", "checkIfAlreadyClassified");
		exit(1);
	} 
}

/**************************************************************************************************/
