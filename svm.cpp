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


LabelPair make_label_pair(const Label& one, const Label& two) {
    LabelVector labelPair;
    labelPair.push_back(one);
    labelPair.push_back(two);
    return labelPair;
}

// the discriminant member function returns +1 or -1
int SVM::discriminant(const Observation& observation) {
    // d is the discriminant function
    double d = b;
    for ( int i = 0; i < y.size(); i++ ) {
        d += y[i]*a[i]*std::inner_product(observation.begin(), observation.end(), x[i].second->begin(), 0.0);
    }
    return d > 0.0 ? 1 : -1;
}

Label SVM::classify(const Observation& observation) {
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


MultiClassSVM::MultiClassSVM(const std::vector<SVM*> s) : twoClassSvmList(s.begin(), s.end()) {
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
Label MultiClassSVM::classify(const Observation& observation) {
    std::map<Label, int> votes;
    for ( int i = 0; i < twoClassSvmList.size(); i++ ) {
        Label predictedLabel = twoClassSvmList[i]->classify(observation);
        votes[predictedLabel]++;
    }
    std::vector<std::pair<Label, int> > tally(votes.begin(), votes.end());
    std::sort(tally.begin(), tally.end(), fewerVotes);
    //for ( int i = 0; i < tally.size(); i++) {
    //    std::cout << tally[i].first << " has " << tally[i].second << " votes" << std::endl;
    //}
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

class MaxIterationsExceeded : public std::exception {
    virtual const char* what() const throw() {
        return "maximum iterations exceeded during SMO";
    }
} maxIterationsExceeded;

//  The train method implements Sequential Minimal Optimization as described in
//  "Support Vector Machine Solvers" by Bottou and Lin.
//
//  SmoTrainer::train releases a pointer to an SVM into the wild so we must be
//  careful about handling the LabeledObservationVector....  Must create a copy
//  of those labeled vectors???
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
    LinearKernelFunction kernelFunction;
    // this is inefficient in general
    std::vector<std::vector<double> > K(observationCount, std::vector<double>(observationCount));
    for ( int u = 0; u < observationCount; u++ ) {
        for ( int v = 0; v < observationCount; v++ ) {
            K[u][v] = kernelFunction.similarity(
                *twoClassLabeledObservationVector[u].second,
                *twoClassLabeledObservationVector[v].second
            );
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
            // what happens if we just go with what we've got instead of throwing an exception?
            // things work pretty well for the most part
            // might be better to look at lambda???
            std::cout << "iteration limit reached" << std::endl;
            break;
            //throw maxIterationsExceeded;
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
    if (verbose) std::cout << "find w" << std::endl;
    std::vector<double> w(twoClassLabeledObservationVector[0].second->size(), 0.0);
    double b = 0.0;
    for ( int i = 0; i < y.size(); i++ ) {
        if (verbose) std::cout << "alpha[" << i << "] = " << a[i] << std::endl;
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

    // be careful about passing twoClassLabeledObservationVector - what if this vector
    // is deleted???
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
        for (LabelSet::const_iterator i = labelSet.begin(); i != labelSet.end(); i++) {
            std::cout << "    label " << *i << std::endl;
        }
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
OneVsOneMultiClassSvmTrainer::OneVsOneMultiClassSvmTrainer(const LabeledObservationVector& lo) :
    labeledObservations(lo) {

    buildLabelSet(labelSet, labeledObservations);
    buildLabelToLabeledObservationVector(labelToLabeledObservationVector, labeledObservations);
    buildLabelPairSet(labelPairSet, labeledObservations);
    standardizeObservations(labeledObservations);
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
        LabelPair labelPair(2);
        labelPair[0] = label;
        for (LabelVector::const_iterator i = labelStack.begin(); i != labelStack.end(); i++) {
            labelPair[1] = *i;
            labelPairSet.insert(
                //std::make_pair(label, *i)
                labelPair
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
        //std::cout << "mean of feature " << feature << " is " << mean << std::endl;
        //std::cout << "std of feature " << feature << " is " << standardDeviation << std::endl;
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
    LabelVector labelList(labelSet.begin(), labelSet.end());
    int deK = 3;
    KFoldLabeledObservationsDivider kFoldDevEvalDivider(deK, labeledObservations);
    kFoldDevEvalDivider.start(labelList);
    const LabeledObservationVector& developmentObservations = kFoldDevEvalDivider.getTrainingData();
    const LabeledObservationVector& evaluationObservations  = kFoldDevEvalDivider.getTestingData();

    // determine hyperparameters by cross validation
    // fill a map with parameter lists
    ParameterRangeMap hyperparameterMap;
    ParameterRange cList;
    // TODO: populate this list from an argument
    cList.push_back(0.01);
    cList.push_back(1.0);
    cList.push_back(10.0);
    hyperparameterMap.insert(std::make_pair("smo_c", cList));

    std::vector<SVM*> twoClassSvmList;
    SmoTrainer smoTrainer;
    LabelPairSet::iterator labelPair;
    for (labelPair = labelPairSet.begin(); labelPair != labelPairSet.end(); labelPair++) {
        // generate training and testing data for this label pair
        Label label0 = (*labelPair)[0];
        Label label1 = (*labelPair)[1];
        std::cout << "training SVM on labels " << label0 << " and " << label1 << std::endl;
        std::cout << "    label pair size is " << labelPair->size() << std::endl;

        // loop on kernel functions and kernel function parameters
        ParameterSetBuilder p(hyperparameterMap);

        // K should be an argument
        int K = 5;
        KFoldLabeledObservationsDivider kFoldLabeledObservationsDivider(K, developmentObservations);
        double bestC = -1.0; // sentinel
        double bestScore = 0.0;
        double meanScoreOverKFolds = 0.0;
        for ( ParameterMapVector::const_iterator hp = p.getParameterSetList().begin(); hp != p.getParameterSetList().end(); hp++ ) {
            // we need to calculate the mean score over the k-folds for the current C
            // oh how I love to calculate mean on-line
            double online_mean_n = 0.0;
            double online_mean_score = 0.0;
            meanScoreOverKFolds = -1.0;  // means we failed to train a SVM

            smoTrainer.setParameters(*hp);

            for ( kFoldLabeledObservationsDivider.start(*labelPair); !kFoldLabeledObservationsDivider.end(); kFoldLabeledObservationsDivider.next() ) {
                const LabeledObservationVector& kthTwoClassTrainingFold = kFoldLabeledObservationsDivider.getTrainingData();
                const LabeledObservationVector& kthTwoClassTestingFold = kFoldLabeledObservationsDivider.getTestingData();
                std::cout << "fold " << kFoldLabeledObservationsDivider.getFoldNumber() << " training data has " << kthTwoClassTrainingFold.size() << " labeled observations" << std::endl;
                std::cout << "fold " << kFoldLabeledObservationsDivider.getFoldNumber() << " testing data has " << kthTwoClassTestingFold.size() << " labeled observations" << std::endl;
                try {
                    SVM* evaluationSvm = smoTrainer.train(kthTwoClassTrainingFold);
                    double score = evaluationSvm->score(kthTwoClassTestingFold);
                    std::cout << "score on fold " << kFoldLabeledObservationsDivider.getFoldNumber() << " test data with C = "<< smoTrainer.getC() << " : " << score << std::endl;
                    online_mean_n += 1.0;
                    double online_mean_delta = score - online_mean_score;
                    online_mean_score += online_mean_delta / online_mean_n;
                    meanScoreOverKFolds = online_mean_score;

                    delete evaluationSvm;
                }
                catch ( std::exception& e ) {
                    std::cout << "exception: " << e.what() << std::endl;
                    std::cout << "    on fold " << kFoldLabeledObservationsDivider.getFoldNumber() << " failed to train SVM with C = " << smoTrainer.getC() << std::endl;
                }
            }
            std::cout << "done with cross validation on C = " << smoTrainer.getC() << std::endl;
            std::cout << "    mean score over " << kFoldLabeledObservationsDivider.getFoldNumber() << " folds is " << meanScoreOverKFolds << std::endl;
            if ( meanScoreOverKFolds > bestScore ) {
                bestC = smoTrainer.getC();
                bestScore = meanScoreOverKFolds;
            }
            else if ( meanScoreOverKFolds < 0.0 ) {
                std::cout << "failed to train SVM with C = " << smoTrainer.getC() << std::endl;
            }
        }
        std::cout << "done with cross validation on all C values" << std::endl;
        if ( bestC < 0.0 ) {
            std::cout << "failed to train SVM on labels " << label0 << " and " << label1 << std::endl;
            throw std::exception();
        }
        else {
            std::cout << "trained SVM on labels " << label0 << " and " << label1 << std::endl;
            std::cout << "    best C is " << bestC << std::endl;

            LabeledObservationVector twoClassDevelopmentObservations;
            std::remove_copy_if(
                developmentObservations.begin(),
                developmentObservations.end(),
                std::back_inserter(twoClassDevelopmentObservations),
                [&](const LabeledObservation& o){
                    return !((o.first == label0) || (o.first == label1));
                }
            );
            std::cout << "training final SVM with C = " << bestC << std::endl;
            std::cout << "training final SVM with " << twoClassDevelopmentObservations.size() << " labeled observations" << std::endl;
            smoTrainer.setC(bestC);
            SVM* svm = smoTrainer.train(twoClassDevelopmentObservations);
            twoClassSvmList.push_back(svm);
        }
    }
    MultiClassSVM* mc = new MultiClassSVM(twoClassSvmList);
    double score = mc->score(evaluationObservations);
    std::cout << "multiclass SVM score: " << score << std::endl;

    return mc;
}
