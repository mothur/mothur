//
//  rftreenode.cpp
//  Mothur
//
//  Created by Sarah Westcott on 10/2/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "rftreenode.hpp"
#include "utils.hpp"

/***********************************************************************/
RFTreeNode::RFTreeNode(vector< vector<int> > bootstrappedTrainingSamples,
                       vector<int> globalDiscardedFeatureIndices,
                       int numFeatures,
                       int numSamples,
                       int numOutputClasses,
                       int generation,
                       int nodeId,
                       float featureStandardDeviationThreshold)

            : bootstrappedTrainingSamples(bootstrappedTrainingSamples),
            globalDiscardedFeatureIndices(globalDiscardedFeatureIndices),
            numFeatures(numFeatures),
            numSamples(numSamples),
            numOutputClasses(numOutputClasses),
            generation(generation),
            isLeaf(false),
            outputClass(-1),
            nodeId(nodeId),
            testSampleMisclassificationCount(0),
            splitFeatureIndex(-1),
            splitFeatureValue(-1),
            splitFeatureEntropy(-1.0),
            ownEntropy(-1.0),
            featureStandardDeviationThreshold(featureStandardDeviationThreshold),
            bootstrappedFeatureVectors(numFeatures, vector<int>(numSamples, 0)),
            bootstrappedOutputVector(numSamples, 0),
            leftChildNode(NULL),
            rightChildNode(NULL),
            parentNode(NULL) {
                
    m = MothurOut::getInstance();
    
    for (int i = 0; i < numSamples; i++) {    // just doing a simple transpose of the matrix
        if (m->getControl_pressed()) { break; }
        for (int j = 0; j < numFeatures; j++) { bootstrappedFeatureVectors[j][i] = bootstrappedTrainingSamples[i][j]; }
    }
    
    for (int i = 0; i < numSamples; i++) { if (m->getControl_pressed()) { break; }
        bootstrappedOutputVector[i] = bootstrappedTrainingSamples[i][numFeatures]; }
    
    createLocalDiscardedFeatureList();
    updateNodeEntropy();
}
/***********************************************************************/
int RFTreeNode::createLocalDiscardedFeatureList(){
    try {
        Utils util;
        for (int i = 0; i < numFeatures; i++) {
                // TODO: need to check if bootstrappedFeatureVectors == numFeatures, in python code we are using bootstrappedFeatureVectors instead of numFeatures
            if (m->getControl_pressed()) { return 0; } 
            vector<int>::iterator it = find(globalDiscardedFeatureIndices.begin(), globalDiscardedFeatureIndices.end(), i);
            if (it == globalDiscardedFeatureIndices.end()) {                           // NOT FOUND
                double standardDeviation = util.getStandardDeviation(bootstrappedFeatureVectors[i]);  
                if (standardDeviation <= featureStandardDeviationThreshold) { localDiscardedFeatureIndices.push_back(i); }
            }
        }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "RFTreeNode", "createLocalDiscardedFeatureList");
        exit(1);
    }  
}
/***********************************************************************/
int RFTreeNode::updateNodeEntropy() {
    try {
        
        vector<int> classCounts(numOutputClasses, 0);
        for (int i = 0; i < bootstrappedOutputVector.size(); i++) {
            classCounts[bootstrappedOutputVector[i]]++;
        }
        int totalClassCounts = accumulate(classCounts.begin(), classCounts.end(), 0);
        double nodeEntropy = 0.0;
        for (int i = 0; i < classCounts.size(); i++) {
            if (m->getControl_pressed()) { return 0; }
            if (classCounts[i] == 0) continue;
            double probability = (double)classCounts[i] / (double)totalClassCounts;
            nodeEntropy += -(probability * log2(probability));
        }
        ownEntropy = nodeEntropy;
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "RFTreeNode", "updateNodeEntropy");
        exit(1);
    } 
}

/***********************************************************************/
