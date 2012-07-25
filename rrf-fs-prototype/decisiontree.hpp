  //
  //  decisiontree.hpp
  //  rrf-fs-prototype
  //
  //  Created by Abu Zaher Faridee on 5/28/12.
  //  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
  //

#ifndef rrf_fs_prototype_decisiontree_hpp
#define rrf_fs_prototype_decisiontree_hpp

#include <vector>
#include <map>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <limits>
#include "macros.h"
#include "treenode.hpp"
#include "abstractdecisiontree.hpp"

using namespace std;

struct VariableRankDescendingSorter {
  bool operator() (vector<int> first, vector<int> second){ return first[1] > second[1]; }
};

class DecisionTree: public AbstractDecisionTree{
  
public:
  
  DecisionTree(vector< vector<int> > baseDataSet,
               vector<int> globalDiscardedFeatureIndices,
               OptimumFeatureSubsetSelector optimumFeatureSubsetSelector,
               string treeSplitCriterion)
  
  : AbstractDecisionTree(baseDataSet,
                         globalDiscardedFeatureIndices,
                         optimumFeatureSubsetSelector,
                         treeSplitCriterion),
  
  variableImportanceList(numFeatures, 0){
    
      // TODO self.outOfBagEstimates = {}
    
//		self.buildDecisionTree()
    
    createBootStrappedSamples();
    buildDecisionTree();
  }
  
  ~DecisionTree(){
    if (rootNode != NULL){ delete rootNode; }
  }
  
  void calcTreeVariableImportanceAndError() {
#ifdef DEBUG_MODE
    DEBUGMSG_LOCATION;
#endif
    int numCorrect;
    double treeErrorRate;
    calcTreeErrorRate(numCorrect, treeErrorRate);
    
    PRINT_VAR(bootstrappedTestSamples.size());
    PRINT_VAR(numCorrect);
    PRINT_VAR(treeErrorRate);
    
    for (unsigned i = 0; i < numFeatures; i++) {
      vector< vector<int> > randomlySampledTestData = randomlyShuffleAttribute(bootstrappedTestSamples, i);
      int numCorrectAfterShuffle = 0;
      for (unsigned j = 0; j < randomlySampledTestData.size(); j++) {
        vector<int> shuffledSample = randomlySampledTestData[j];
        int actualSampleOutputClass = shuffledSample[numFeatures];
        int predictedSampleOutputClass = evaluateSample(shuffledSample);
        if (actualSampleOutputClass == predictedSampleOutputClass) { numCorrectAfterShuffle++; }
        variableImportanceList[i] += (numCorrect - numCorrectAfterShuffle);
      }
    }
    
      // TODO: do we need to save the variableRanks in the DecisionTree, do we need it later?
    vector< vector<int> > variableRanks;
    for (unsigned i = 0; i < variableImportanceList.size(); i++) {
      if (variableImportanceList[i] > 0) {
          // TODO: is there a way to optimize the follow line's code?
        vector<int> variableRank(2, 0);
        variableRank[0] = i; variableRank[1] = variableImportanceList[i];
      }
    }
    VariableRankDescendingSorter variableRankDescendingSorter;
    sort(variableRanks.begin(), variableRanks.end(), variableRankDescendingSorter);
    PRINT_VAR(variableRanks);
  }
  
  int evaluateSample(vector<int> testSample) {
    TreeNode *node = rootNode;
    while (true) {
      if (node->checkIsLeaf() == true) { return node->getOutputClass(); }
      int sampleSplitFeatureValue = testSample[node->getSplitFeatureIndex()];
      if (sampleSplitFeatureValue < node->getSplitFeatureValue()) { node = node->getLeftChildNode(); }
      else { node = node->getRightChildNode(); } 
    }
    return 0;
  }
  
  void calcTreeErrorRate(int& numCorrect, double& treeErrorRate){
    numCorrect = 0;
    for (unsigned i = 0; i < bootstrappedTestSamples.size(); i++) {
      
      vector<int> testSample = bootstrappedTestSamples[i];
      int testSampleIndex = bootstrappedTestSampleIndices[i];
      
      int actualSampleOutputClass = testSample[numFeatures];
      int predictedSampleOutputClass = evaluateSample(testSample);
      
      if (actualSampleOutputClass == predictedSampleOutputClass) { numCorrect++; } 
      
      outOfBagEstimates[testSampleIndex] = predictedSampleOutputClass;
    }
    
    treeErrorRate = 1 - ((double)numCorrect / (double)bootstrappedTestSamples.size());    
  }
  
    // TODO: optimize the algo, instead of transposing two time, we can extarct the feature,
    // shuffle it and then re-insert in the original place, thus iproving runnting time
  vector< vector<int> > randomlyShuffleAttribute(vector< vector<int> > samples, int featureIndex) {
    
    vector< vector<int> > featureVectors(numFeatures, vector<int>(samples.size(), 0));
    for (unsigned i = 0; i < samples.size(); i++) { for (unsigned j = 0; j < numFeatures; i++) { featureVectors[j][i] =  samples[i][j]; } }
    
    random_shuffle(featureVectors[featureIndex].begin(), featureVectors[featureIndex].end());
    
    vector< vector<int> > shuffledSample(samples.size(), vector<int>(numFeatures, 0));
    for (unsigned i = 0; i < shuffledSample.size(); i++) { for (unsigned j = 0; j < numFeatures; i++) { shuffledSample[i][j] = featureVectors[j][i]; } }
    
    return shuffledSample;
  }
  
protected:
  
private:
  
  void buildDecisionTree(){
    int generation = 0;
    rootNode = new TreeNode(bootstrappedTrainingSamples, globalDiscardedFeatureIndices, numFeatures, numSamples, numOutputClasses, generation);
    splitRecursively(rootNode);
#ifdef DEBUG_MODE
    printTree(rootNode, "root");
#endif
  }
  
  void splitRecursively(TreeNode* rootNode){
    
#ifdef DEBUG_MODE
    DEBUGMSG_LOCATION;
#endif
    
    if (rootNode->getNumSamples() < 2){      
#ifdef DEBUG_MODE
      DEBUGMSG("Already classified: Case 1");
#endif      
      rootNode->setIsLeaf(true);
      rootNode->setOutputClass(rootNode->getBootstrappedTrainingSamples()[0][rootNode->getNumFeatures()]);
      return;
    }
    
    int classifiedOutputClass;
    bool isAlreadyClassified = checkIfAlreadyClassified(rootNode, classifiedOutputClass);    
    if (isAlreadyClassified == true){
#ifdef DEBUG_MODE
      DEBUGMSG("Already classified: Case 2");
#endif
      rootNode->setIsLeaf(true);
      rootNode->setOutputClass(classifiedOutputClass);
      return;
    }
    
    vector<int> featureSubsetIndices = selectFeatureSubsetRandomly(globalDiscardedFeatureIndices, rootNode->getLocalDiscardedFeatureIndices());
    rootNode->setFeatureSubsetIndices(featureSubsetIndices);
    
#ifdef DEBUG_MODE
    DEBUGMSG_VAR(globalDiscardedFeatureIndices);
    DEBUGMSG_VAR(featureSubsetIndices);
#endif
    
    findAndUpdateBestFeatureToSplitOn(rootNode);
    
    vector< vector<int> > leftChildSamples;
    vector< vector<int> > rightChildSamples;
    getSplitPopulation(rootNode, leftChildSamples, rightChildSamples);
    
      // TODO: need to write code to clear this memory
    TreeNode* leftChildNode = new TreeNode(leftChildSamples, globalDiscardedFeatureIndices, numFeatures, leftChildSamples.size(), numOutputClasses, rootNode->getGeneration() + 1);
    TreeNode* rightChildNode = new TreeNode(rightChildSamples, globalDiscardedFeatureIndices, numFeatures, rightChildSamples.size(), numOutputClasses, rootNode->getGeneration() + 1);
    
    rootNode->setLeftChildNode(leftChildNode);
    leftChildNode->setParentNode(rootNode);
    
    rootNode->setRightChildNode(rightChildNode);
    rightChildNode->setParentNode(rootNode);
    
      // TODO: This recursive split can be parrallelized later
    splitRecursively(leftChildNode);
    splitRecursively(rightChildNode);
  }
  
    // TODO: finish implementation of findAndUpdateBestFeatureToSplitOn()
  void findAndUpdateBestFeatureToSplitOn(TreeNode* node){
    
#ifdef DEBUG_MODE
    DEBUGMSG_LOCATION;
#endif
    
    vector< vector<int> > bootstrappedFeatureVectors = node->getBootstrappedFeatureVectors();
    vector<int> bootstrappedOutputVector = node->getBootstrappedOutputVector();
    vector<int> featureSubsetIndices = node->getFeatureSubsetIndices();
    
    vector<double> featureSubsetEntropies;
    vector<int> featureSubsetSplitValues;
    vector<double> featureSubsetIntrinsicValues;
    vector<double> featureSubsetGainRatios;
    
    for (unsigned i = 0; i < featureSubsetIndices.size(); i++) {
      int tryIndex = featureSubsetIndices[i];
      
#ifdef DEBUG_MODE
      cout << "trying feature of index:" << tryIndex << endl;
#endif
      double featureMinEntropy;
      int featureSplitValue;
      double featureIntrinsicValue;
      
      getMinEntropyOfFeature(bootstrappedFeatureVectors[tryIndex], bootstrappedOutputVector, featureMinEntropy, featureSplitValue, featureIntrinsicValue);
      
      featureSubsetEntropies.push_back(featureMinEntropy);
      featureSubsetSplitValues.push_back(featureSplitValue);
      featureSubsetIntrinsicValues.push_back(featureIntrinsicValue);
      
      double featureInformationGain = node->getOwnEntropy() - featureMinEntropy;
      double featureGainRatio = (double)featureInformationGain / (double)featureIntrinsicValue;
      featureSubsetGainRatios.push_back(featureGainRatio);
      
    }
    
#ifdef DEBUG_MODE
    DEBUGMSG_VAR(featureSubsetEntropies);
    DEBUGMSG_VAR(featureSubsetSplitValues);
    DEBUGMSG_VAR(featureSubsetIntrinsicValues);
    DEBUGMSG_VAR(featureSubsetGainRatios);
#endif

    vector<double>::iterator minEntropyIterator = min_element(featureSubsetEntropies.begin(), featureSubsetEntropies.end());
    vector<double>::iterator maxGainRatioIterator = max_element(featureSubsetGainRatios.begin(), featureSubsetGainRatios.end());
    double featureMinEntropy = *minEntropyIterator;
    double featureMaxGainRatio = *maxGainRatioIterator;

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
    
#ifdef DEBUG_MODE
    DEBUGMSG_VAR(bestFeatureToSplitOnIndex);
    DEBUGMSG_VAR(bestFeatureSplitValue);
    DEBUGMSG_VAR(bestFeatureSplitEntropy);
    if (treeSplitCriterion == "gainRatio") { DEBUGMSG_VAR(featureMaxGainRatio); }
#endif
    
    node->setSplitFeatureIndex(featureSubsetIndices[bestFeatureToSplitOnIndex]);
    node->setSplitFeatureValue(bestFeatureSplitValue);
    node->setSplitFeatureEntropy(bestFeatureSplitEntropy);
  }
  
  vector<int> selectFeatureSubsetRandomly(vector<int> globalDiscardedFeatureIndices, vector<int> localDiscardedFeatureIndices){
#ifdef DEBUG_MODE
    DEBUGMSG_LOCATION;
#endif
    vector<int> featureSubsetIndices;
    
    vector<int> combinedDiscardedFeatureIndices;
    combinedDiscardedFeatureIndices.insert(combinedDiscardedFeatureIndices.end(), globalDiscardedFeatureIndices.begin(), globalDiscardedFeatureIndices.end());
    combinedDiscardedFeatureIndices.insert(combinedDiscardedFeatureIndices.end(), localDiscardedFeatureIndices.begin(), localDiscardedFeatureIndices.end());
    
    sort(combinedDiscardedFeatureIndices.begin(), combinedDiscardedFeatureIndices.end());
    
    int numberOfRemainingSuitableFeatures = numFeatures - combinedDiscardedFeatureIndices.size();
    int currentFeatureSubsetSize = numberOfRemainingSuitableFeatures < optimumFeatureSubsetSize ? numberOfRemainingSuitableFeatures : optimumFeatureSubsetSize;
    
    while (featureSubsetIndices.size() < currentFeatureSubsetSize) {
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
    
    return featureSubsetIndices;
  }
  
    // TODO: printTree() needs a check if correct
  void printTree(TreeNode* treeNode, string caption){
    
    string tabs = "";
    for (unsigned i = 0; i < treeNode->getGeneration(); i++) { tabs += "\t"; }
    
    if (treeNode != NULL && treeNode->checkIsLeaf() == false){
      cout << tabs << caption << " [ gen: " << treeNode->getGeneration() << " ] ( " << treeNode->getSplitFeatureValue() << " < X" << treeNode->getSplitFeatureIndex() << " )" << endl;
      printTree(treeNode->getLeftChildNode(), "leftChild");
      printTree(treeNode->getRightChildNode(), "rightChild");
    }else {
      cout << tabs << caption << " [ gen: " << treeNode->getGeneration() + " ] ( classified to: " << treeNode->getOutputClass() << ", samples: " << treeNode->getNumSamples() << " )";
    }
    
  }
  
  void purgeDataSetsFromTree() {
    purgeTreeNodesDataRecursively(rootNode);
  }
  
  void purgeTreeNodesDataRecursively(TreeNode* treeNode) {
    treeNode->bootstrappedTrainingSamples.clear();
    treeNode->bootstrappedFeatureVectors.clear();
    treeNode->bootstrappedOutputVector.clear();
    treeNode->localDiscardedFeatureIndices.clear();
    treeNode->globalDiscardedFeatureIndices.clear();
    
    if (treeNode->leftChildNode != NULL) { purgeTreeNodesDataRecursively(treeNode->leftChildNode); }
    if (treeNode->rightChildNode != NULL) { purgeTreeNodesDataRecursively(treeNode->rightChildNode); }
  }
  
  vector<int> variableImportanceList;
};

#endif
