//
//  svm.cpp
//  support vector machine
//
//  Created by Joshua Lynch on 6/19/2013.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//
#include <iostream>
#include <numeric>
#include <stack>
#include <utility>

#include "svm.hpp"


MultiClassSVM::MultiClassSVM() {
}

MultiClassSVM::~MultiClassSVM() {
}

SmoTrainer::SmoTrainer() : C(1.0){
}

SmoTrainer::~SmoTrainer() {
}

//  The train method implements Sequential Minimal Optimization as described in
//  "Support Vector Machine Solvers" by Bottou and Lin.
//
//  Here is a helpful example showing the details of this implementation.
//  Consider a dataset of 4 observations, each having 2 features:
//    observation  feature 1  feature 2  label    y (class)
//      1            0.0        0.0        blue     -1.0
//      2            1.0        0.0        green    +1.0
//      3            0.0        1.0        blue     -1.0
//      4            1.0        1.0        green    +1.0
//
//  The dual coefficients a and the gradient g are initially:
//    a = (0.0. 0.0, 0.0, 0.0)
//    g = (1.0, 1.0, 1.0, 1.0)
//  The box constraint vectors A and B are:
//    A = (-1.0, 0.0, -1.0, 0.0)
//    B = ( 0.0, 1.0,  0.0, 1.0)
//  The (linear) kernel matrix:
//      1 2 3 4
//    1 0 0 0 0
//    2 0 4 0 4
//    3 0 0 1 1
//    4 0 4 1 5
//
//  First time through the loop:
//        A    B    y    a    g    ya   yg   ya<B yg[n]>yg[i] i
//    1  -1.0  0.0 -1.0  0.0  1.0  0.0 -1.0   f        f      0
//    2   0.0  1.0 +1.0  0.0  1.0  0.0  1.0   t        t      1
//    3  -1.0  0.0 -1.0  0.0  1.0  0.0 -1.0   f        f      1
//    4   0.0  1.0 +1.0  0.0  1.0  0.0  1.0   t        f      1

SVM* SmoTrainer::train(const ObservationVector& twoClassObservationVector, const LabelVector& twoClassLabelVector) {
    const int observationCount = twoClassObservationVector.size();
    const int featureCount = twoClassObservationVector[0]->size();
    std::cout << "observation count : " << observationCount << std::endl;
    std::cout << "feature count     : " << featureCount << std::endl;
    // dual coefficients
    std::vector<double> a(observationCount, 0.0);
    // gradient
    std::vector<double> g(observationCount, 1.0);
    // convert the labels to -1.0,+1.0
    std::vector<double> y(observationCount);
    std::cout << "assign numeric labels" << std::endl;
    assignNumericLabels(y, twoClassLabelVector);
    std::cout << "assign A and B" << std::endl;
    std::vector<double> A(observationCount);
    std::vector<double> B(observationCount);
    for ( int n = 0; n < observationCount; n++ ) {
        if ( y[n] == +1.0) {
            A[n] = 0.0;
            B[n] = C;
        }
        else {
            A[n] = -C;
            B[n] = 0;
        }
        std::cout << n << " " << A[n] << " " << B[n] << std::endl;
    }
    std::cout << "assign K" << std::endl;
    // this is inefficient in general
    std::vector<std::vector<double> > K(observationCount, std::vector<double>(observationCount));
    for ( int u = 0; u < observationCount; u++ ) {
        for ( int v = 0; v < observationCount; v++ ) {
            K[u][v] = std::inner_product(
                twoClassObservationVector[u]->begin(),
                twoClassObservationVector[u]->end(),
                twoClassObservationVector[v]->begin(),
                0.0);
            std::cout << "u = " << u << " v = " << v << " K = " << K[u][v] << std::endl;
        }
    }
    std::cout << "train" << std::endl;
    std::vector<double> u(3);
    std::vector<double> ya(observationCount);
    std::vector<double> yg(observationCount);
    while ( true ) {
        int i = 0;
        int j = 0;
        for ( int n = 0; n < observationCount; n++ ) {
            ya[n] = y[n] * a[n];
            yg[n] = y[n] * g[n];
            if ( ya[n] < B[n] && yg[n] > yg[i] ) {
                i = n;
            }
            if ( A[n] < ya[n] && yg[n] < yg[j] ) {
                j = n;
            }
            std::cout << "n = " << n << std::endl;
            std::cout << "i = " << i << " yg[i] = " << yg[i] << std::endl;
            std::cout << "j = " << j << " yg[j] = " << yg[j] << std::endl;
        }
        // maximum violating pair is i,j
        //std::cout << "stopping criterion:" << std::endl;

        if ( yg[i] <= yg[j] ) {
            break;
        }
        u[0] = B[i] - ya[i];
        u[1] = ya[j] - A[j];
        u[2] = (yg[i] - yg[j]) / (K[i][i]+K[j][j]-2.0*K[i][j]);
        double lambda = *std::min_element(u.begin(), u.end());
        std::cout << "lambda: " << lambda << std::endl;
        for ( int k = 0; k < featureCount; k++ ) {
            g[k] -= lambda * y[k] * K[i][k] + lambda * y[k] * K[j][k];
        }
        a[i] += y[i] * lambda;
        a[j] -= y[j] * lambda;
    }

    return new SVM(y, a);
}

// For SVM we need to assign numeric labels of -1.0 and +1.0.
// This method populates the y vector argument with -1.0 and +1.0
// corresponding to the two classes in the labelVector argument.
void SmoTrainer::assignNumericLabels(std::vector<double>& y, const LabelVector& labelVector) {
    // it would be nice if we assign -1.0 and +1.0 consistently for each pair of labels
	// I think the set will always be traversed in sorted order so we should get this for free
    LabelSet labelSet(labelVector.begin(), labelVector.end());
    LabelVector uniqueLabels(labelSet.begin(), labelSet.end());
    if (labelSet.size() != 2) {
        // throw an exception
        throw std::exception();
    }
    else {
        // should make certain y has correct size?
        // this loop will populate the y vector with
        // -1.0 for uniqueLabel[0] and +1.0 for uniqueLabel[1]
        for ( int i = 0; i < labelVector.size(); i++ ) {
            y[i] = (labelVector[i] == uniqueLabels[0] ? -1.0 : 1.0);
        }
    }
    std::cout << "returning from assign numeric labels" << std::endl;
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
