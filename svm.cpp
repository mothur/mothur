//
//  svm.cpp
//  support vector machine
//
//  Created by Joshua Lynch on 6/19/2013.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//
#include <algorithm>
#include <iostream>
#include <limits>
#include <numeric>
#include <stack>
#include <utility>

#include "svm.hpp"

// the classify method should return a class label
Label classifier::empty("empty");

// the discriminant member function returns +1 or -1
int SVM::discriminant(const FeatureVector& observation) {
    // d is the discriminant function
    double d = b;
    for ( int i = 0; i < y.size(); i++ ) {
        d += y[i]*a[i]*std::inner_product(observation.begin(), observation.end(), x[i].second->begin(), 0.0);
    }
    //std::cout << "d = " << d << std::endl;
    return d > 0.0 ? 1 : -1;
}

Label SVM::classify(const FeatureVector& observation) {
    int d = discriminant(observation);
    return discriminantToLabel[d];
}

double SVM::score(const LabeledObservationVector& twoClassLabeledObservationVector) {
    double s = 0.0;
    for ( int i = 0; i < twoClassLabeledObservationVector.size(); i++ ) {
        Label predicted_label = classify(*twoClassLabeledObservationVector[i].second);
        if ( predicted_label == twoClassLabeledObservationVector[i].first ) {
            s = s + 1.0;
        }
        else {

        }
    }
    return s / double(twoClassLabeledObservationVector.size());
}


MultiClassSVM::MultiClassSVM(const std::vector<SVM*> s) : twoClassSvmList(s) {
}

MultiClassSVM::~MultiClassSVM() {
    for ( int i = 0; i < twoClassSvmList.size(); i++ ) {
        delete twoClassSvmList[i];
    }
}

bool fewerVotes(const std::pair<Label, int>& p, const std::pair<Label, int>& q) {
    return p.second < q.second;
}

// TODO: handle ties
Label MultiClassSVM::classify(const FeatureVector& observation) {
    std::map<Label, int> votes;
    for ( int i = 0; i < twoClassSvmList.size(); i++ ) {
        Label predictedLabel = twoClassSvmList[i]->classify(observation);
        votes[predictedLabel]++;
    }
    std::vector<std::pair<Label, int> > tally(votes.begin(), votes.end());
    std::sort(tally.begin(), tally.end(), fewerVotes);
    for ( int i = 0; i < tally.size(); i++) {
        //std::cout << tally[i].first << " has " << tally[i].second << " votes" << std::endl;
    }
    std::pair<Label, int> winner = *tally.rbegin();
    Label winningLabel = winner.first;
    int winningVoteCount = winner.second;
    //std::cout << "winner is " << winningLabel << " with " << winningVoteCount << " votes" << std::endl;
    return winningLabel;
}

double MultiClassSVM::score(const LabeledObservationVector& twoClassLabeledObservationVector) {
    double s = 0.0;
    for ( int i = 0; i < twoClassLabeledObservationVector.size(); i++ ) {
        Label predicted_label = classify(*twoClassLabeledObservationVector[i].second);
        if ( predicted_label == twoClassLabeledObservationVector[i].first ) {
            s = s + 1.0;
        }
        else {

        }
    }
    return s / double(twoClassLabeledObservationVector.size());
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
//      2            2.0        0.0        green    +1.0
//      3            0.0        1.0        blue     -1.0
//      4            2.0        1.0        green    +1.0
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
//    n   A    B    y    a    g    ya   yg   ya<B  i  yg[n]>yg[i] A<ya yg[n]<yg[j] i
//    1  -1.0  0.0 -1.0  0.0  1.0  0.0 -1.0    f   0       f       t        f 0
//    2   0.0  1.0 +1.0  0.0  1.0  0.0  1.0    t       t       f   1
//    3  -1.0  0.0 -1.0  0.0  1.0  0.0 -1.0    f       f       t   1
//    4   0.0  1.0 +1.0  0.0  1.0  0.0  1.0    t       f       f   1

SVM* SmoTrainer::train(const LabeledObservationVector& twoClassLabeledObservationVector) {
    const int observationCount = twoClassLabeledObservationVector.size();
    const int featureCount = twoClassLabeledObservationVector[0].second->size();
    bool verbose = false;
    if (verbose) std::cout << "observation count : " << observationCount << std::endl;
    if (verbose) std::cout << "feature count     : " << featureCount << std::endl;
    // dual coefficients
    std::vector<double> a(observationCount, 0.0);
    // gradient
    std::vector<double> g(observationCount, 1.0);
    // convert the labels to -1.0,+1.0
    std::vector<double> y(observationCount);
    if (verbose) std::cout << "assign numeric labels" << std::endl;
    NumericClassToLabel discriminantToLabel;
    assignNumericLabels(y, twoClassLabeledObservationVector, discriminantToLabel);
    if (verbose) std::cout << "assign A and B" << std::endl;
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
        if (verbose) std::cout << n << " " << A[n] << " " << B[n] << std::endl;
    }
    if (verbose) std::cout << "assign K" << std::endl;
    // this is inefficient in general
    std::vector<std::vector<double> > K(observationCount, std::vector<double>(observationCount));
    for ( int u = 0; u < observationCount; u++ ) {
        for ( int v = 0; v < observationCount; v++ ) {
            K[u][v] = std::inner_product(
                twoClassLabeledObservationVector[u].second->begin(),
                twoClassLabeledObservationVector[u].second->end(),
                twoClassLabeledObservationVector[v].second->begin(),
                0.0);
            //std::cout << "u = " << u << " v = " << v << " K = " << K[u][v] << std::endl;
            if (verbose) std::cout << " " << K[u][v];
        }
        if (verbose) std::cout << std::endl;
    }
    int m = 0;
    std::cout << "train" << std::endl;
    std::vector<double> u(3);
    std::vector<double> ya(observationCount);
    std::vector<double> yg(observationCount);
    while ( true ) {
        m++;
        int i = 0;
        int j = 0;
        double yg_max = std::numeric_limits<double>::min();
        double yg_min = std::numeric_limits<double>::max();
        if (verbose) std::cout << "m = " << m << std::endl;
        for ( int k = 0; k < observationCount; k++ ) {
            ya[k] = y[k] * a[k];
            yg[k] = y[k] * g[k];
        }
        if (verbose) std::cout << "yg =";
        for ( int k = 0; k < observationCount; k++ ) {
            //std::cout << A[k] << " " << B[k] << " " << y[k] << " " << a[k] << " " << g[k] << " " << ya[k] << " " << yg[k] << std::endl;
            if (verbose) std::cout << " " << yg[k];
        }
        if (verbose) std::cout << std::endl;

        for ( int k = 0; k < observationCount; k++ ) {
            //ya[k] = y[k] * a[k];
            //yg[k] = y[k] * g[k];
            if ( ya[k] < B[k] && yg[k] > yg_max ) {
                yg_max = yg[k];
                i = k;
            }
            if ( A[k] < ya[k] && yg[k] < yg_min ) {
                yg_min = yg[k];
                j = k;
            }
            //std::cout << "n = " << n << std::endl;
            //std::cout << ya[k] << " " << yg[k] << std::endl;
            //std::cout << "j = " << j << " yg[j] = " << yg[j] << std::endl;
        }
        // maximum violating pair is i,j
        if (verbose) std::cout << "maximal violating pair:" << std::endl;
        if (verbose) std::cout << "  i = " << i << " ";
        for ( int feature = 0; feature < featureCount; feature++ ) {
            if (verbose) std::cout << twoClassLabeledObservationVector[i].second->at(feature) << " ";
        };
        if (verbose) std::cout << std::endl;
        if (verbose) std::cout << "  j = " << j << " ";
        for ( int feature = 0; feature < featureCount; feature++ ) {
            if (verbose) std::cout << twoClassLabeledObservationVector[j].second->at(feature) << " ";
        };
        if (verbose) std::cout << std::endl;
        //std::cout << "stopping criterion:" << std::endl;

        // parameterize this
        if ( m > 1000) {
            throw std::exception();
        }

        if ( yg[i] <= yg[j] ) {
            break;
        }
        u[0] = B[i] - ya[i];
        u[1] = ya[j] - A[j];
        //std::cout << "(" << yg[i] << " - " << yg[j] << ") / (" << K[i][i] << " + " << K[j][j] << " - 2.0 * " << K[i][j] << ")" << std::endl;
        u[2] = (yg[i] - yg[j]) / (K[i][i]+K[j][j]-2.0*K[i][j]);
        if (verbose) std::cout << "directions: (" << u[0] << "," << u[1] << "," << u[2] << ")" << std::endl;
        double lambda = *std::min_element(u.begin(), u.end());
        if (verbose) std::cout << "lambda: " << lambda << std::endl;
        for ( int k = 0; k < observationCount; k++ ) {
            g[k] += (-lambda * y[k] * K[i][k] + lambda * y[k] * K[j][k]);
        }
        a[i] += y[i] * lambda;
        a[j] -= y[j] * lambda;
    }


    // at this point the optimal a's have been found
    // now use them to find w and b
    std::vector<double> w(twoClassLabeledObservationVector[0].second->size(), 0.0);
    double b = 0.0;
    for ( int i = 0; i < y.size(); i++ ) {
        //std::cout << "alpha[" << j << "] = " << a[j] << std::endl;
        for ( int j = 0; j < w.size(); j++ ) {
            w[j] += a[i] * y[i] * twoClassLabeledObservationVector[i].second->at(j);
        }
        if ( A[i] < a[i] && a[i] < B[i] ) {
            b = yg[i];
            if (verbose) std::cout << "b = " << b << std::endl;
        }
    }

    for ( int i = 0; i < w.size(); i++ ) {
        if (verbose) std::cout << "w[" << i << "] = " << w[i] << std::endl;
    }

    return new SVM(y, a, twoClassLabeledObservationVector, b, discriminantToLabel);
}

// For SVM we need to assign numeric labels of -1.0 and +1.0.
// This method populates the y vector argument with -1.0 and +1.0
// corresponding to the two classes in the labelVector argument.
void SmoTrainer::assignNumericLabels(std::vector<double>& y, const LabeledObservationVector& labeledObservationVector, NumericClassToLabel& discriminantToLabel) {
    // it would be nice if we assign -1.0 and +1.0 consistently for each pair of labels
	// I think the set will always be traversed in sorted order so we should get this for free
    LabelSet labelSet;
    OneVsOneMultiClassSvmTrainer::buildLabelSet(labelSet, labeledObservationVector);
    LabelVector uniqueLabels(labelSet.begin(), labelSet.end());
    if (labelSet.size() != 2) {
        // throw an exception
        std::cout << "unexpected label set size " << labelSet.size() << std::endl;
        throw std::exception();
    }
    else {
        // should make certain y has correct size?
        // this loop will populate the y vector with
        // -1.0 for uniqueLabel[0] and +1.0 for uniqueLabel[1]
        for ( int i = 0; i < labeledObservationVector.size(); i++ ) {
            y[i] = (labeledObservationVector[i].first == uniqueLabels[0] ? -1.0 : 1.0);
        }
        discriminantToLabel[-1] = uniqueLabels[0];
        discriminantToLabel[+1] = uniqueLabels[1];
    }
}


//
// OneVsOneMultiClassSvmTrainer
//
// An instance of OneVsOneMultiClassSvmTrainer is intended to work with a single set of data
// to produce a single instance of MultiClassSVM.  That's why observations and labels go in to
// the constructor.
//OneVsOneMultiClassSvmTrainer::OneVsOneMultiClassSvmTrainer(const ObservationVector& o, const LabelVector& l) :
//observations(o), observationLabels(l) {
OneVsOneMultiClassSvmTrainer::OneVsOneMultiClassSvmTrainer(const LabeledObservationVector& lo) :
    labeledObservations(lo) {

    buildLabelSet(labelSet, labeledObservations);
    buildLabelToLabeledObservationVector(labelToLabeledObservationVector, labeledObservations);
    buildLabelPairSet(labelPairSet, labeledObservations);
    standardizeObservations(labeledObservations);
}

OneVsOneMultiClassSvmTrainer::~OneVsOneMultiClassSvmTrainer() {
    std::cout << "~OneVsOneMultiClassSvmTrainer" << std::endl;
    // no need to delete anything
}

void OneVsOneMultiClassSvmTrainer::buildLabelSet(LabelSet& labelSet, const LabeledObservationVector& labeledObservationVector) {
    for (LabeledObservationVector::const_iterator i = labeledObservationVector.begin(); i != labeledObservationVector.end(); i++) {
        labelSet.insert(i->first);
    }
}

void OneVsOneMultiClassSvmTrainer::buildLabelToLabeledObservationVector(LabelToLabeledObservationVector& labelToLabeledObservationVector, const LabeledObservationVector& labeledObservationVector) {
    std::cout << "buildLabelToLabelObservationVector" << std::endl;
    for ( LabeledObservationVector::const_iterator j = labeledObservationVector.begin(); j != labeledObservationVector.end(); j++ ) {
        labelToLabeledObservationVector[j->first].push_back(*j);
    }
}

//  This function uses the LabeledObservationVector argument to populate the LabelPairSet
//  argument with pairs of labels.  For example, if labeledObservationVector looks like this:
//    [ ("blue", x), ("green", y), ("red", z) ]
//  then the labelPairSet will be populated with the following label pairs:
//    ("blue", "green"), ("blue", "red"), ("green", "red")
//  The order of labels in the pairs is determined by the ordering of labels in the temporary
//  LabelSet.  By default this order will be ascending.  However, labels are taken off the
//  temporary labelStack in reverse order, so the labelStack is initialized with reverse iterators.
//  In the end our label pairs will be in sorted order.
void OneVsOneMultiClassSvmTrainer::buildLabelPairSet(LabelPairSet& labelPairSet, const LabeledObservationVector& labeledObservationVector) {
    //std::cout << "buildLabelPairSet" << std::endl;
    LabelSet labelSet;
    buildLabelSet(labelSet, labeledObservationVector);
    LabelVector labelStack(labelSet.rbegin(), labelSet.rend());
    while (labelStack.size() > 1) {
        Label label = labelStack.back();
        labelStack.pop_back();
        for (LabelVector::const_iterator i = labelStack.begin(); i != labelStack.end(); i++) {
            labelPairSet.insert(
                std::make_pair(label, *i)
            );
        }
    }
}

void OneVsOneMultiClassSvmTrainer::standardizeObservations(const LabeledObservationVector& observations) {
    std::cout << "standardizeObservations" << std::endl;
    // online method for mean and variance
    for ( FeatureVector::size_type feature = 0; feature < observations[0].second->size(); feature++ ) {
        double n = 0.0;
        double mean = 0.0;
        double M2 = 0.0;
        for ( ObservationVector::size_type observation = 0; observation < observations.size(); observation++ ) {
            n += 1.0;
            double x = observations[observation].second->at(feature);
            double delta = x - mean;
            mean += delta / n;
            M2 += delta * (x - mean);
        }
        double variance = M2 / (n - 1.0);
        double standardDeviation = sqrt(variance);
        std::cout << "mean of feature " << feature << " is " << mean << std::endl;
        std::cout << "std of feature " << feature << " is " << standardDeviation << std::endl;
        // normalize the feature
        for ( ObservationVector::size_type observation = 0; observation < observations.size(); observation++ ) {
            observations[observation].second->at(feature) = (observations[observation].second->at(feature) - mean ) / standardDeviation;
        }
    }
}

// this method has too many arguments
// split the observations and labels into a training set and a testing set
// for now 2/3 and 1/3
void OneVsOneMultiClassSvmTrainer::appendTrainingAndTestingData(
        Label label,
        const LabeledObservationVector& observationsWithLabel,
        LabeledObservationVector& trainingObservations,
        LabeledObservationVector& testingObservations) {

    for (int i = 0; i < observationsWithLabel.size(); i++) {
        if (i % 3 == 0) {
            testingObservations.push_back(observationsWithLabel[i]);
        }
        else {
            trainingObservations.push_back(observationsWithLabel[i]);
        }
    }
}

MultiClassSVM* OneVsOneMultiClassSvmTrainer::train() {
    // divide the data into development and evaluation sets
    LabelToLabeledObservationVector labeledDevelopmentObservations;
    //LabeledObservations labeledEvaluationObservations;
    LabeledObservationVector developmentObservations;
    LabeledObservationVector evaluationObservations;
    //LabelVector developmentLabels;
    //LabelVector evaluationLabels;
    for (LabelSet::iterator label = labelSet.begin(); label != labelSet.end(); label++) {
        appendTrainingAndTestingData(*label, labelToLabeledObservationVector[*label], developmentObservations, evaluationObservations);
    }
    buildLabelToLabeledObservationVector(labeledDevelopmentObservations, developmentObservations);

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

    std::vector<SVM*> twoClassSvmList;
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
    for (labelPair = labelPairSet.begin(); labelPair != labelPairSet.end(); labelPair++) {
        // generate training and testing data for this label pair
        Label label_0 = labelPair->first;
        Label label_1 = labelPair->second;
        std::cout << "training SVM on labels " << label_0 << " and " << label_1 << std::endl;
        LabeledObservationVector twoClassTrainingVector;

        LabeledObservationVector twoClassTestingVector;

        appendTrainingAndTestingData(label_0, labeledDevelopmentObservations[label_0], twoClassTrainingVector, twoClassTestingVector);
        appendTrainingAndTestingData(label_1, labeledDevelopmentObservations[label_1], twoClassTrainingVector, twoClassTestingVector);

        SVM* bestSvm = NULL;
        double bestScore = 0.0;
        for (hp = cList.begin(); hp != cList.end(); ++hp) {
            double C = *hp;
            smoTrainer.setC(C);
            try {
                SVM* evaluationSvm = smoTrainer.train(twoClassTrainingVector);
                double score = evaluationSvm->score(twoClassTestingVector);
                std::cout << "score on test data for C = "<< C << " : " << score << std::endl;
                if ( score > bestScore ) {
                    bestSvm = evaluationSvm;
                }
                else {
                    delete evaluationSvm;
                }
            }
            catch (std::exception& e) {
                std::cout << e.what() << std::endl;
                std::cout << "failed to train SVM with C = " << C << std::endl;
            }
        }
        if ( bestSvm == NULL ) {
            std::cout << "failed to train SVM on labels " << label_0 << " and " << label_1 << std::endl;
            throw std::exception();
        }
        else {
            twoClassSvmList.push_back(bestSvm);
        }
    }

    MultiClassSVM* mc = new MultiClassSVM(twoClassSvmList);
    double score = mc->score(evaluationObservations);
    std::cout << "multiclass SVM score: " << score << std::endl;

    return mc;
}
