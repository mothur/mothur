//
//  svm.cpp
//  support vector machine
//
//  Created by Joshua Lynch on 6/19/2013.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//
#include <iostream>
#include <utility>
#include <stack>

#include "svm.hpp"

SVM::SVM() {
	// TODO Auto-generated constructor stub

}

SVM::~SVM() {
	// TODO Auto-generated destructor stub
}

MultiClassSVM::MultiClassSVM() {
}

MultiClassSVM::~MultiClassSVM() {
}

SmoTrainer::SmoTrainer() : C(0.0){
}

SmoTrainer::~SmoTrainer() {
}

//
// OneVsOneMultiClassSvmTrainer
//
// An instance of OneVsOneMultiClassSvmTrainer is intended to work with a single set of data
// to produce a single instance of MultiClassSVM.  That's why observations and labels go in to
// the constructor.
// the ObservationVector use here is not quite right
OneVsOneMultiClassSvmTrainer::OneVsOneMultiClassSvmTrainer(const ObservationVector& o, const LabelVector& l) :
    observations(o), observationLabels(l) {
	labelSet.clear();
	labelSet.insert(observationLabels.begin(), observationLabels.end());
    // should maybe require a vector of label-observation pairs?
    // insert label, observation vector pairs
    LabelSet::iterator i;
    for (i = labelSet.begin(); i != labelSet.end(); i++) {
        std::cout << "creating observation vector for label " << *i << std::endl;
        //labeledObservations.insert(make_pair(*i, ObservationVector()));
    }

}

OneVsOneMultiClassSvmTrainer::~OneVsOneMultiClassSvmTrainer() {
    // destroy the ObservationVectors
    LabeledObservations::iterator i;
    for (i = labeledObservations.begin(); i != labeledObservations.end(); i++) {
        std::cout << "deleting observation vector for label " << i->first << std::endl;
        //delete i->second;
    }
}

SVM* OneVsOneMultiClassSvmTrainer::train(
        const ObservationVector& observations,
        const LabelVector& observationLabels) {

    // TODO: enumerate the labels, figure out how many two-class classifiers we need
    // insert labeled observations in a map of labels to observation vectors
    LabeledObservations labeledObservations;

    LabelSet labelSet;
    //getLabelSet(labelSet, observationLabels);
    labelSet.clear();
    labelSet.insert(observationLabels.begin(), observationLabels.end());
    // fix this method name it stinks
    //getLabeledObservations(labeledObservations, observations, observationLabels);
    // insert label, observation vector pairs
    //LabelSet::iterator i;
    //for (i = labelSet.begin(); i != labelSet.end(); i++) {
    //    Label label = *i;
    //    labeledObservations.insert(make_pair(label, new ObservationVector()));
    //}
    std::set<LabelPair> labelPairs;
    getLabelPairSet(labelPairs, labelSet);

    // split the observations and labels into a training set and a testing set
    // for now 2/3 and 1/3
    int numObservations = observations.size();
    ObservationVector trainingData;
    LabelVector trainingLabels;
    ObservationVector testingData;
    LabelVector testingLabels;
    // TODO: use pointers or references
    // TODO: be certain the labels are distributed evenly between training and test sets
    for (int i = 0; i < observations.size(); i++) {
        if (i % 3 == 0) {
            testingData.push_back(observations[i]);
            testingLabels.push_back(observationLabels[i]);
        }
        else {
            trainingData.push_back(observations[i]);
            testingLabels.push_back(observationLabels[i]);
        }
    }

    // determine hyperparameters by cross validation
    // fill a map with parameter lists
    std::map<std::string, std::vector<double> > hyperparameterMap;
    typedef std::vector<double> HyperparameterList;
    HyperparameterList cList;
    // TODO: populate this list from an argument
    cList.push_back(0.01);
    cList.push_back(1.0);
    cList.push_back(10.0);
    hyperparameterMap.insert(std::make_pair("C", cList));
    hyperparameterMap.insert(std::make_pair("D", cList)); // fake

    SmoTrainer smoTrainer;
    //std::map<std::string, std::vector<double> >::iterator pos;
    // keep track of the performance of each hyperparameter set with a map of
    //    keys   : map of strings to doubles (hyperparameter names to values)
    //    values : double (classification performance - fraction correct)
    // first build the keys and insert 0.0 values
    // punt
    // first loop on label pairs, second loop on hyperparameters
    LabelPairSet::iterator labelPair;
    HyperparameterList::iterator hp;
    for (labelPair = labelPairs.begin(); labelPair != labelPairs.end(); labelPair++) {
        // generate training and testing data for this label pair
        ObservationVector twoClassTrainingVector;
        LabelVector twoClassTrainingLabelVector;
        ObservationVector twoClassTestingVector;
        LabelVector twoClassTestingLabelVector;
        for (hp = cList.begin(); hp != cList.end(); ++hp) {
            double C = *hp;
            smoTrainer.setC(C);
            SVM* evaluationSvm = smoTrainer.train(twoClassTrainingVector, twoClassTrainingLabelVector);
            double score = evaluationSvm->score(twoClassTestingVector, twoClassTestingLabelVector);
        }
    }
    return NULL;
}

//void OneVsOneMultiClassSvmTrainer::getLabelSet(LabelSet& labelSet, const LabelVector& labelVector) {
//    labelSet.clear();
//    labelSet.insert(labelVector.begin(), labelVector.end());
//}

void OneVsOneMultiClassSvmTrainer::getLabeledObservations(LabeledObservations& labeledObservations, const ObservationVector& observationVector, const LabelVector& labelVector) {
    LabelSet labelSet;
    //getLabelSet(labelSet, labelVector);
    labelSet.clear();
    labelSet.insert(labelVector.begin(), labelVector.end());

    // the LabeledObservation map must deallocate
    labeledObservations.clear();
    LabelSet::iterator i;
    for (i = labelSet.begin(); i != labelSet.end(); i++) {
        ObservationVector* v = new ObservationVector();
        //labeledObservations.insert(make_pair(*i, v));
    }
}

void OneVsOneMultiClassSvmTrainer::getLabelPairSet(LabelPairSet& labelPairSet, const LabelSet& labelSet) {
    labelPairSet.clear();
    LabelVector labelStack(labelSet.begin(), labelSet.end());
    while (labelStack.size() > 1) {
        Label label = labelStack.back();
        labelStack.pop_back();
        for (int i = 0; i < labelStack.size(); i++) {
            labelPairSet.insert(
                std::make_pair(label, labelStack[i])
            );
        }
    }
}
