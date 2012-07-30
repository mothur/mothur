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
  
  friend class RandomForest;
  
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
        
    createBootStrappedSamples();
    buildDecisionTree();
    
  }
  
  virtual ~DecisionTree(){
    deleteTreeNodesRecursively(rootNode);
  }
  
  void calcTreeVariableImportanceAndError() {
    
#ifdef DEBUG_LEVEL_2
    DEBUGMSG_FUNC;
#endif
    
    int numCorrect;
    double treeErrorRate;
    calcTreeErrorRate(numCorrect, treeErrorRate);
    
    PRINT_VAR(bootstrappedTestSamples.size());
    PRINT_VAR(numCorrect);
    PRINT_VAR(treeErrorRate);
        
    for (unsigned i = 0; i < numFeatures; i++) {
        // NOTE: only shuffle the features, never shuffle the output vector
        // so i = 0 and i will be alwaays <= (numFeatures - 1) as the index at numFeatures will denote
        // the feature vector
      vector< vector<int> > randomlySampledTestData = randomlyShuffleAttribute(bootstrappedTestSamples, i);
      
      int numCorrectAfterShuffle = 0;
      for (unsigned j = 0; j < randomlySampledTestData.size(); j++) {
        vector<int> shuffledSample = randomlySampledTestData[j];
        int actualSampleOutputClass = shuffledSample[numFeatures];
        int predictedSampleOutputClass = evaluateSample(shuffledSample);
        if (actualSampleOutputClass == predictedSampleOutputClass) { numCorrectAfterShuffle++; }
      }
      variableImportanceList[i] += (numCorrect - numCorrectAfterShuffle);
    }
    
      // TODO: do we need to save the variableRanks in the DecisionTree, do we need it later?
    vector< vector<int> > variableRanks;
    for (unsigned i = 0; i < variableImportanceList.size(); i++) {
      if (variableImportanceList[i] > 0) {
          // TODO: is there a way to optimize the follow line's code?
        vector<int> variableRank(2, 0);
        variableRank[0] = i; variableRank[1] = variableImportanceList[i];
        variableRanks.push_back(variableRank);
      }
    }
    VariableRankDescendingSorter variableRankDescendingSorter;
    sort(variableRanks.begin(), variableRanks.end(), variableRankDescendingSorter);
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
    
      // NOTE: we need (numFeatures + 1) featureVecotors, the last extra vector is actually outputVector
    vector< vector<int> > featureVectors(samples[0].size(), vector<int>(samples.size(), 0));
        
    for (unsigned i = 0; i < samples.size(); i++) {
      for (unsigned j = 0; j < samples[0].size(); j++) {
        featureVectors[j][i] =  samples[i][j];
      }
    }
    
    random_shuffle(featureVectors[featureIndex].begin(), featureVectors[featureIndex].end());
    
    vector< vector<int> > shuffledSample(samples.size(), vector<int>(samples[0].size(), 0));

    for (unsigned i = 0; i < shuffledSample.size(); i++) {
      for (unsigned j = 0; j < numFeatures; j++) {
        shuffledSample[i][j] = featureVectors[j][i];
      }
    }
    
    return shuffledSample;
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
    
#ifdef DEBUG_LEVEL_2
    DEBUGMSG_FUNC;
#endif
    
    if (rootNode->getNumSamples() < 2){
      
#ifdef DEBUG_LEVEL_2
      DEBUGMSG("Already classified: Case 1");
#endif
      
      rootNode->setIsLeaf(true);
      rootNode->setOutputClass(rootNode->getBootstrappedTrainingSamples()[0][rootNode->getNumFeatures()]);
      return;
    }
    
    int classifiedOutputClass;
    bool isAlreadyClassified = checkIfAlreadyClassified(rootNode, classifiedOutputClass);    
    if (isAlreadyClassified == true){
      
#ifdef DEBUG_LEVEL_2
      DEBUGMSG("Already classified: Case 2");
#endif
      
      rootNode->setIsLeaf(true);
      rootNode->setOutputClass(classifiedOutputClass);
      return;
    }
    
    vector<int> featureSubsetIndices = selectFeatureSubsetRandomly(globalDiscardedFeatureIndices, rootNode->getLocalDiscardedFeatureIndices());
    rootNode->setFeatureSubsetIndices(featureSubsetIndices);
    
#ifdef DEBUG_LEVEL_2
    PRINT_VAR(globalDiscardedFeatureIndices);
    PRINT_VAR(featureSubsetIndices);
#endif
        
    findAndUpdateBestFeatureToSplitOn(rootNode);
    
    vector< vector<int> > leftChildSamples;
    vector< vector<int> > rightChildSamples;
    getSplitPopulation(rootNode, leftChildSamples, rightChildSamples);
    
      // TODO: need to write code to clear this memory
    TreeNode* leftChildNode = new TreeNode(leftChildSamples, globalDiscardedFeatureIndices, numFeatures, (int)leftChildSamples.size(), numOutputClasses, rootNode->getGeneration() + 1);
    TreeNode* rightChildNode = new TreeNode(rightChildSamples, globalDiscardedFeatureIndices, numFeatures, (int)rightChildSamples.size(), numOutputClasses, rootNode->getGeneration() + 1);
    
    rootNode->setLeftChildNode(leftChildNode);
    leftChildNode->setParentNode(rootNode);
    
    rootNode->setRightChildNode(rightChildNode);
    rightChildNode->setParentNode(rootNode);
    
      // TODO: This recursive split can be parrallelized later
    splitRecursively(leftChildNode);
    splitRecursively(rightChildNode);
  }
  
  void findAndUpdateBestFeatureToSplitOn(TreeNode* node){
    
#ifdef DEBUG_LEVEL_2
    DEBUGMSG_FUNC;
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
      
#ifdef DEBUG_LEVEL_4
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
    
#ifdef DEBUG_LEVEL_2
    PRINT_VAR(featureSubsetIndices);
    PRINT_VAR(featureSubsetEntropies);
    PRINT_VAR(featureSubsetSplitValues);
    PRINT_VAR(featureSubsetIntrinsicValues);
    PRINT_VAR(featureSubsetGainRatios);
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
    
#ifdef DEBUG_LEVEL_3
    PRINT_VAR(bestFeatureToSplitOnIndex);
    PRINT_VAR(featureSubsetIndices[bestFeatureToSplitOnIndex]);
    PRINT_VAR(bestFeatureSplitValue);
    PRINT_VAR(bestFeatureSplitEntropy);
    if (treeSplitCriterion == "gainRatio") { PRINT_VAR(featureMaxGainRatio); }
#endif
    
    node->setSplitFeatureIndex(featureSubsetIndices[bestFeatureToSplitOnIndex]);
    node->setSplitFeatureValue(bestFeatureSplitValue);
    node->setSplitFeatureEntropy(bestFeatureSplitEntropy);
  }
  
  vector<int> selectFeatureSubsetRandomly(vector<int> globalDiscardedFeatureIndices, vector<int> localDiscardedFeatureIndices){
    
#ifdef DEBUG_LEVEL_2
    DEBUGMSG_FUNC;
#endif
    
    vector<int> featureSubsetIndices;
    
    vector<int> combinedDiscardedFeatureIndices;
    combinedDiscardedFeatureIndices.insert(combinedDiscardedFeatureIndices.end(), globalDiscardedFeatureIndices.begin(), globalDiscardedFeatureIndices.end());
    combinedDiscardedFeatureIndices.insert(combinedDiscardedFeatureIndices.end(), localDiscardedFeatureIndices.begin(), localDiscardedFeatureIndices.end());
        
    sort(combinedDiscardedFeatureIndices.begin(), combinedDiscardedFeatureIndices.end());
    
#ifdef DEBUG_LEVEL_2
    PRINT_VAR(combinedDiscardedFeatureIndices);
#endif
    
    int numberOfRemainingSuitableFeatures = (int)(numFeatures - combinedDiscardedFeatureIndices.size());
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
    
//#ifdef DEBUG_LEVEL_2
//    PRINT_VAR(featureSubsetIndices);
//#endif
    
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
  
  void deleteTreeNodesRecursively(TreeNode* treeNode) {
    if (treeNode == NULL) { return; }
    deleteTreeNodesRecursively(treeNode->leftChildNode);
    deleteTreeNodesRecursively(treeNode->rightChildNode);
    delete treeNode;
  }
  
  vector<int> variableImportanceList;
  map<int, int> outOfBagEstimates;
};

#endif
