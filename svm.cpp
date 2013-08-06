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

#define RANGE(X) X, X + sizeof(X)/sizeof(double)

const std::string LinearKernelFunction::MapKey     = "LinearKernel";
const std::string LinearKernelFunction::MapKey_Constant = "LinearKernel_Constant";
const double defaultLinearConstantRangeArray[5] = {-10.0, -1.0, 0.0, 1.0, 10.0};
const ParameterRange LinearKernelFunction::defaultConstantRange = ParameterRange(RANGE(defaultLinearConstantRangeArray));

const std::string RbfKernelFunction::MapKey        = "RbfKernel";
const std::string RbfKernelFunction::MapKey_Gamma  = "RbfKernel_Gamma";
const double defaultRbfGammaRangeArray[] = {0.0001, 0.001, 0.01, 0.1, 1.0, 10.0, 100.0};
const ParameterRange RbfKernelFunction::defaultGammaRange = ParameterRange(RANGE(defaultRbfGammaRangeArray));

const std::string PolynomialKernelFunction::MapKey          = "PolynomialKernel";
const std::string PolynomialKernelFunction::MapKey_Constant = "PolynomialKernel_Constant";
const std::string PolynomialKernelFunction::MapKey_Degree   = "PolynomialKernel_Degree";

const double defaultPolynomialConstantRangeArray[] = {-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0};
const ParameterRange PolynomialKernelFunction::defaultConstantRange = ParameterRange(RANGE(defaultPolynomialConstantRangeArray));
const double defaultPolynomialDegreeRangeArray[] = {2.0, 3.0, 4.0};
const ParameterRange PolynomialKernelFunction::defaultDegreeRange = ParameterRange(RANGE(defaultPolynomialDegreeRangeArray));

const std::string SigmoidKernelFunction::MapKey          = "SigmoidKernel";
const std::string SigmoidKernelFunction::MapKey_Alpha    = "SigmoidKernel_Alpha";
const std::string SigmoidKernelFunction::MapKey_Constant = "SigmoidKernel_Constant";

const double defaultSigmoidAlphaRangeArray[] = {1.0, 2.0};
const ParameterRange SigmoidKernelFunction::defaultAlphaRange = ParameterRange(RANGE(defaultSigmoidAlphaRangeArray));
const double defaultSigmoidConstantRangeArray[] = {1.0, 2.0};
const ParameterRange SigmoidKernelFunction::defaultConstantRange = ParameterRange(RANGE(defaultSigmoidConstantRangeArray));

const std::string SmoTrainer::MapKey_C = "SmoTrainer_C";
const double defaultSmoTrainerCRangeArray[] = {0.0001, 0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0};
const ParameterRange SmoTrainer::defaultCRange = ParameterRange(RANGE(defaultSmoTrainerCRangeArray));


LabelPair buildLabelPair(const Label& one, const Label& two) {
    LabelVector labelPair(2);
    labelPair[0] = one;
    labelPair[1] = two;
    return labelPair;
}

// Dividing a dataset into training and testing sets while maintaing equal
// representation of all classes is done using a LabelToLabeledObservationVector.
// This container is used to divide datasets into groups of LabeledObservations
// having the same label.  For example, given a LabeledObservationVector like
//     ["blue",  [1.0, 2.0, 3.0]]
//     ["green", [3.0, 4.0, 5.0]]
//     ["blue",  [2,0, 3.0. 4.0]]
//     ["green", [4.0, 5.0, 6.0]]
// the corresponding LabelToLabeledObservationVector looks like
//     "blue",  [["blue",  [1.0, 2.0, 3.0]], ["blue",  [2,0, 3.0. 4.0]]]
//     "green", [["green", [3.0, 4.0, 5.0]], ["green", [4.0, 5.0, 6.0]]]
void buildLabelToLabeledObservationVector(LabelToLabeledObservationVector& labelToLabeledObservationVector, const LabeledObservationVector& labeledObservationVector) {
    for ( LabeledObservationVector::const_iterator j = labeledObservationVector.begin(); j != labeledObservationVector.end(); j++ ) {
        labelToLabeledObservationVector[j->first].push_back(*j);
    }
}


//
// SVM member functions
//
// the discriminant member function returns +1 or -1
int SVM::discriminant(const Observation& observation) {
    // d is the discriminant function
    double d = b;
    for ( int i = 0; i < y.size(); i++ ) {
        d += y[i]*a[i]*std::inner_product(observation.begin(), observation.end(), x[i].second->begin(), 0.0);
    }
    return d > 0.0 ? 1 : -1;
}

// the score member function classifies each labeled observation from the
// argument and returns the percent of correct classifications
double SVM::score(const LabeledObservationVector& twoClassLabeledObservationVector) {
    double s = 0.0;
    for (LabeledObservationVector::const_iterator i = twoClassLabeledObservationVector.begin(); i != twoClassLabeledObservationVector.end(); i++) {
        Label predicted_label = classify(*(i->second));
        if ( predicted_label == i->first ) {
            s = s + 1.0;
        }
        else {

        }
    }
    return s / double(twoClassLabeledObservationVector.size());
}

// The fewerVotes function is used to find the maximum vote
// tally in MultiClassSVM::classify.  This function returns true
// if the first element (number of votes for the first label) is
// less than the second element (number of votes for the second label).
bool fewerVotes(const std::pair<Label, int>& p, const std::pair<Label, int>& q) {
    return p.second < q.second;
}


class MultiClassSvmClassificationTie : public std::exception {
public:
    MultiClassSvmClassificationTie(LabelVector& t, int c) : tiedLabels(t), tiedVoteCount(c) {}
    ~MultiClassSvmClassificationTie() throw() {}

    virtual const char* what() const throw() {
        return "classification tie";
    }

private:
    const LabelVector tiedLabels;
    const int tiedVoteCount;
};

Label MultiClassSVM::classify(const Observation& observation) {
    std::map<Label, int> labelToVoteCount;
    for ( int i = 0; i < twoClassSvmList.size(); i++ ) {
        Label predictedLabel = twoClassSvmList[i]->classify(observation);
        labelToVoteCount[predictedLabel]++;
    }
    std::pair<Label, int> winner = *max_element(labelToVoteCount.begin(), labelToVoteCount.end(), fewerVotes);
    LabelVector winningLabels;
    winningLabels.push_back(winner.first);
    for ( std::map<Label, int>::const_iterator i = labelToVoteCount.begin(); i != labelToVoteCount.end(); i++ ) {
        if ( i->second == winner.second && i->first != winner.first ) {
            winningLabels.push_back(i->first);
        }
    }
    if ( winningLabels.size() == 1) {
        // we have a winner
    }
    else {
        // we have a tie
        throw MultiClassSvmClassificationTie(winningLabels, winner.second);
    }

    return winner.first;
}

double MultiClassSVM::score(const LabeledObservationVector& multiClassLabeledObservationVector) {
    double s = 0.0;
    for (LabeledObservationVector::const_iterator i = multiClassLabeledObservationVector.begin(); i != multiClassLabeledObservationVector.end(); i++) {
        //std::cout << "classifying observation with label " << i->first << std::endl;
        try {
            Label predicted_label = classify(*(i->second));
            if ( predicted_label == i->first ) {
                s = s + 1.0;
            }
            else {

            }
        }
        catch ( MultiClassSvmClassificationTie& m ) {
            std::cout << "classification tie for observation with label " << i->first << std::endl;
        }
    }
    return s / double(multiClassLabeledObservationVector.size());
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
SVM* SmoTrainer::train(KernelFunction* kernelFunction, const LabeledObservationVector& twoClassLabeledObservationVector) {
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
    std::vector<std::vector<double> > K(observationCount, std::vector<double>(observationCount, std::numeric_limits<double>::quiet_NaN()));
    /*
    for ( int u = 0; u < observationCount; u++ ) {
        for ( int v = 0; v < observationCount; v++ ) {
            K[u][v] = kernelFunction->similarity(
                *twoClassLabeledObservationVector[u].second,
                *twoClassLabeledObservationVector[v].second
            );
            if (verbose) std::cout << " " << K[u][v];
        }
        if (verbose) std::cout << std::endl;
    }
    */
    int m = 0;
    std::vector<double> u(3);
    std::vector<double> ya(observationCount);
    std::vector<double> yg(observationCount);
    double lambda = std::numeric_limits<double>::max();
    while ( true ) {
        m++;
        int i = 0; // 0
        int j = 0; // 0
        double yg_max = std::numeric_limits<double>::min();
        double yg_min = std::numeric_limits<double>::max();
        if (verbose) std::cout << "m = " << m << std::endl;
        for ( int k = 0; k < observationCount; k++ ) {
            ya[k] = y[k] * a[k];
            yg[k] = y[k] * g[k];
        }
        if (verbose) {
            std::cout << "yg =";
            for ( int k = 0; k < observationCount; k++ ) {
                //std::cout << A[k] << " " << B[k] << " " << y[k] << " " << a[k] << " " << g[k] << " " << ya[k] << " " << yg[k] << std::endl;
                std::cout << " " << yg[k];
            }
            std::cout << std::endl;
        }

        for ( int k = 0; k < observationCount; k++ ) {
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
        if (verbose) {
            std::cout << "maximal violating pair: " << i << " " << j << std::endl;
            std::cout << "  i = " << i << " features: ";
            for ( int feature = 0; feature < featureCount; feature++ ) {
                std::cout << twoClassLabeledObservationVector[i].second->at(feature) << " ";
            };
            std::cout << std::endl;
            std::cout << "  j = " << j << " features: ";
            for ( int feature = 0; feature < featureCount; feature++ ) {
                std::cout << twoClassLabeledObservationVector[j].second->at(feature) << " ";
            };
            std::cout << std::endl;
        }

        // parameterize this
        if ( m > 1000 ) { //1000
            // what happens if we just go with what we've got instead of throwing an exception?
            // things work pretty well for the most part
            // might be better to look at lambda???
            std::cout << "iteration limit reached with lambda = " << lambda << std::endl;
            break;
        }

        // using lambda to break is a good performance enhancement
        if ( yg[i] <= yg[j] or lambda < 0.0001) {
            break;
        }
        u[0] = B[i] - ya[i];
        u[1] = ya[j] - A[j];
        // calculate similarities for i and j if we have not done so yet
        //std::cout << "K[" << i << "][0] = " << K[i][0] << std::endl;
        if (std::isnan(K[i][0])) {
            //std::cout << "calculating row " << i << " of K" << std::endl;
            //for ( int u = 0; u < observationCount; u++ ) {
                for ( int v = 0; v < observationCount; v++ ) {
                    K[i][v] = kernelFunction->similarity(
                        *twoClassLabeledObservationVector[i].second,
                        *twoClassLabeledObservationVector[v].second);
                    //K[v][i] = K[i][v];
                    //if (verbose) std::cout << " " << K[u][v];
                }
                //if (verbose) std::cout << std::endl;
            //}
            //std::cout << "    done" << std::endl;
        }
        //std::cout << "K[" << j << "][0] = " << K[j][0] << std::endl;
        if (std::isnan(K[j][0])) {
            //std::cout << "calculating row " << j << " of K" << std::endl;
            //for ( int u = 0; u < observationCount; u++ ) {
                for ( int v = 0; v < observationCount; v++ ) {
                    K[j][v] = kernelFunction->similarity(
                        *twoClassLabeledObservationVector[j].second,
                        *twoClassLabeledObservationVector[v].second);
                    //K[v][j] = K[j][v];
                    //if (verbose) std::cout << " " << K[u][v];
                }
                //if (verbose) std::cout << std::endl;
            //}
            //std::cout << "    done" << std::endl;
        }

        u[2] = (yg[i] - yg[j]) / (K[i][i]+K[j][j]-2.0*K[i][j]);
        if (verbose) std::cout << "directions: (" << u[0] << "," << u[1] << "," << u[2] << ")" << std::endl;
        lambda = *std::min_element(u.begin(), u.end());
        if (verbose) std::cout << "lambda: " << lambda << std::endl;
        for ( int k = 0; k < observationCount; k++ ) {
            g[k] += (-lambda * y[k] * K[i][k] + lambda * y[k] * K[j][k]);
        }
        if (verbose) {
            std::cout << "g =";
            for ( int k = 0; k < observationCount; k++ ) {
                std::cout << " " << g[k];
            }
            std::cout << std::endl;
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

    if (verbose) {
        for ( int i = 0; i < w.size(); i++ ) {
            std::cout << "w[" << i << "] = " << w[i] << std::endl;
        }
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
    buildLabelSet(labelSet, labeledObservationVector);
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

void getDefaultKernelParameterRangeMap(KernelParameterRangeMap& kernelParameterRangeMap) {
    ParameterRangeMap linearParameterRangeMap;
    linearParameterRangeMap[SmoTrainer::MapKey_C] = SmoTrainer::defaultCRange;
    linearParameterRangeMap[LinearKernelFunction::MapKey_Constant] = LinearKernelFunction::defaultConstantRange;

    ParameterRangeMap rbfParameterRangeMap;
    rbfParameterRangeMap[SmoTrainer::MapKey_C] = SmoTrainer::defaultCRange;
    rbfParameterRangeMap[RbfKernelFunction::MapKey_Gamma] = RbfKernelFunction::defaultGammaRange;

    ParameterRangeMap polynomialParameterRangeMap;
    polynomialParameterRangeMap[SmoTrainer::MapKey_C] = SmoTrainer::defaultCRange;
    polynomialParameterRangeMap[PolynomialKernelFunction::MapKey_Constant] = PolynomialKernelFunction::defaultConstantRange;
    polynomialParameterRangeMap[PolynomialKernelFunction::MapKey_Degree] = PolynomialKernelFunction::defaultDegreeRange;

    kernelParameterRangeMap[LinearKernelFunction::MapKey] = linearParameterRangeMap;
    kernelParameterRangeMap[RbfKernelFunction::MapKey] = rbfParameterRangeMap;
    kernelParameterRangeMap[PolynomialKernelFunction::MapKey] = polynomialParameterRangeMap;
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

void buildLabelSet(LabelSet& labelSet, const LabeledObservationVector& labeledObservationVector) {
    for (LabeledObservationVector::const_iterator i = labeledObservationVector.begin(); i != labeledObservationVector.end(); i++) {
        labelSet.insert(i->first);
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
    for ( Observation::size_type feature = 0; feature < observations[0].second->size(); feature++ ) {
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

// The LabelMatchesEither functor is used only in a call to remove_copy_if in the
// OneVsOneMultiClassSvmTrainer::train method.  It returns true if the labeled
// observation argument has the same label as either of the two label arguments.
class LabelMatchesEither {
public:
    LabelMatchesEither(const Label& _label0, const Label& _label1) : label0(_label0), label1(_label1) {}

    bool operator() (const LabeledObservation& o) {
        return !((o.first == label0) || (o.first == label1));
    }

private:
    const Label& label0;
    const Label& label1;
};

MultiClassSVM* OneVsOneMultiClassSvmTrainer::train(const KernelParameterRangeMap& kernelParameterRangeMap) {
    // first divide the data into a 'development' set for tuning hyperparameters
    // and an 'evaluation' set for measuring performance
    LabelVector labelList(labelSet.begin(), labelSet.end());
    int deK = 3;
    KFoldLabeledObservationsDivider kFoldDevEvalDivider(deK, labeledObservations);
    kFoldDevEvalDivider.start();
    const LabeledObservationVector& developmentObservations = kFoldDevEvalDivider.getTrainingData();
    const LabeledObservationVector& evaluationObservations  = kFoldDevEvalDivider.getTestingData();

    std::vector<SVM*> twoClassSvmList;
    SmoTrainer smoTrainer;
    LabelPairSet::iterator labelPair;
    for (labelPair = labelPairSet.begin(); labelPair != labelPairSet.end(); labelPair++) {
        // generate training and testing data for this label pair
        Label label0 = (*labelPair)[0];
        Label label1 = (*labelPair)[1];
        //std::cout << "training SVM on labels " << label0 << " and " << label1 << std::endl;
        //std::cout << "    label pair size is " << labelPair->size() << std::endl;

        double bestScore = 0.0;
        ParameterMap bestParameterMap;
        std::string bestKernelFunctionKey;
        // K should be an argument
        int K = 5;
        LabeledObservationVector twoClassDevelopmentObservations;
        LabelMatchesEither labelMatchesEither(label0, label1);
        std::remove_copy_if(
            developmentObservations.begin(),
            developmentObservations.end(),
            std::back_inserter(twoClassDevelopmentObservations),
            labelMatchesEither
            //[&](const LabeledObservation& o){
            //    return !((o.first == label0) || (o.first == label1));
            //}
        );
        KFoldLabeledObservationsDivider kFoldLabeledObservationsDivider(K, twoClassDevelopmentObservations);
        // loop on kernel functions and kernel function parameters
        for ( KernelParameterRangeMap::const_iterator kmap = kernelParameterRangeMap.begin(); kmap != kernelParameterRangeMap.end(); kmap++ ) {
            std::string kernelFunctionKey = kmap->first;
            std::cout << "training with kernel " << kmap->first << std::endl;
            KernelFunction* kernelFunction = KernelFactory::getKernelFunctionForKey(kmap->first);
            ParameterSetBuilder p(kmap->second);
            for (ParameterMapVector::const_iterator hp = p.getParameterSetList().begin(); hp != p.getParameterSetList().end(); hp++) {
                kernelFunction->setParameters(*hp);
                smoTrainer.setParameters(*hp);
                for ( ParameterMap::const_iterator i = hp->begin(); i != hp->end(); i++ ) {
                    std::cout << "    " << i->first << ":" << i->second << std::endl;
                }
                double score = trainOnKFolds(smoTrainer, kernelFunction, kFoldLabeledObservationsDivider);
                if ( score > bestScore ) {
                    bestScore = score;
                    bestParameterMap = *hp;
                    bestKernelFunctionKey = kernelFunctionKey;
                }
            }
            delete kernelFunction;
        }

        std::cout << "done with cross validation on all parameter values" << std::endl;
        if ( bestScore == 0.0 ) {
            std::cout << "failed to train SVM on labels " << label0 << " and " << label1 << std::endl;
            throw std::exception();
        }
        else {
            std::cout << "trained SVM on labels " << label0 << " and " << label1 << std::endl;
            std::cout << "    best score is " << bestScore << std::endl;
            std::cout << "    best parameters are " << std::endl;
            for ( ParameterMap::const_iterator p = bestParameterMap.begin(); p != bestParameterMap.end(); p++ ) {
                std::cout << "        "  << p->first << " : " << p->second << std::endl;
            }

            LabelMatchesEither labelMatchesEither(label0, label1);
            LabeledObservationVector twoClassDevelopmentObservations;
            std::remove_copy_if(
                developmentObservations.begin(),
                developmentObservations.end(),
                std::back_inserter(twoClassDevelopmentObservations),
                labelMatchesEither
                //[&](const LabeledObservation& o){
                //    return !((o.first == label0) || (o.first == label1));
                //}
            );
            //std::cout << "training final SVM with C = " << bestC << std::endl;
            std::cout << "training final SVM with " << twoClassDevelopmentObservations.size() << " labeled observations" << std::endl;
            for ( ParameterMap::const_iterator i = bestParameterMap.begin(); i != bestParameterMap.end(); i++ ) {
                std::cout << "    " << i->first << ":" << i->second << std::endl;
            }
            KernelFunction* kernelFunction = KernelFactory::getKernelFunctionForKey(bestKernelFunctionKey);
            kernelFunction->setParameters(bestParameterMap);
            smoTrainer.setParameters(bestParameterMap);
            SVM* svm = smoTrainer.train(kernelFunction, twoClassDevelopmentObservations);
            std::cout << "done training final SVM" << std::endl;
            twoClassSvmList.push_back(svm);
        }
    }

    std::cout << "building MultiClassSvm" << std::endl;
    MultiClassSVM* mc = new MultiClassSVM(twoClassSvmList);
    std::cout << "done building MultiClassSvm" << std::endl;
    double score = mc->score(evaluationObservations);
    std::cout << "multiclass SVM score: " << score << std::endl;

    return mc;
}

double OneVsOneMultiClassSvmTrainer::trainOnKFolds(SmoTrainer& smoTrainer, KernelFunction* kernelFunction, KFoldLabeledObservationsDivider& kFoldLabeledObservationsDivider) {
    double meanScoreOverKFolds = 0.0;
    double online_mean_n = 0.0;
    double online_mean_score = 0.0;
    meanScoreOverKFolds = -1.0;  // means we failed to train a SVM

    for ( kFoldLabeledObservationsDivider.start(); !kFoldLabeledObservationsDivider.end(); kFoldLabeledObservationsDivider.next() ) {
        const LabeledObservationVector& kthTwoClassTrainingFold = kFoldLabeledObservationsDivider.getTrainingData();
        const LabeledObservationVector& kthTwoClassTestingFold = kFoldLabeledObservationsDivider.getTestingData();
        //std::cout << "fold " << kFoldLabeledObservationsDivider.getFoldNumber() << " training data has " << kthTwoClassTrainingFold.size() << " labeled observations" << std::endl;
        //std::cout << "fold " << kFoldLabeledObservationsDivider.getFoldNumber() << " testing data has " << kthTwoClassTestingFold.size() << " labeled observations" << std::endl;
        try {
            SVM* evaluationSvm = smoTrainer.train(kernelFunction, kthTwoClassTrainingFold);
            double score = evaluationSvm->score(kthTwoClassTestingFold);
            std::cout << "score on fold " << kFoldLabeledObservationsDivider.getFoldNumber() << " of test data is " << score << std::endl;
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
    if ( meanScoreOverKFolds == 0.0 ) {
        std::cout << "failed to train SVM with C = " << smoTrainer.getC() << std::endl;
    }
    return meanScoreOverKFolds;
}

// Only the linear kernel can be used here.
/*
FeatureLabelVector SvmRfe::getOrderedFeatureList(const LabeledObservationVector& labeledObservationVector, const KernelParameterRangeMap& linearKernelParameterRangeMap) {
    OneVsOneMultiClassSvmTrainer t(labeledObservationVector);
    FeatureLabelVector rankedFeatureVector;
    while ( rankedFeatureVector.size() < labeledObservationVector.at(0).second->size() ) {
        MultiClassSVM* s = t.train(kernelParameterRangeMap);
        // calculate the 'ranking criterion' for each feature using each binary svm
        //for each feature i
        //    for each svm r
        //        c_i += (w_r_i^2)
        //delete s;
    }
}
*/
