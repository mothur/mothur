//
//  svm.hpp
//  support vector machine
//
//  Created by Joshua Lynch on 6/19/2013.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#ifndef svm_hpp_
#define svm_hpp_

#include <algorithm>
#include <cmath>
#include <exception>
#include <map>
#include <set>
#include <string>
#include <vector>

// a FeatureVector is an Observation ....
typedef std::vector<double> FeatureVector;
typedef std::vector<double> Observation;
typedef std::vector<Observation*> ObservationVector;
typedef std::string Label;
typedef std::vector<Label> LabelVector;
typedef std::set<std::string> LabelSet;

//typedef std::pair<Label,Label> LabelPair;
typedef std::vector<Label> LabelPair;

typedef std::set<LabelPair> LabelPairSet;
typedef std::pair<Label, Observation*> LabeledObservation;
typedef std::vector<LabeledObservation> LabeledObservationVector;
typedef std::map<Label, LabeledObservationVector> LabelToLabeledObservationVector;
typedef std::map<int, Label> NumericClassToLabel;

LabelPair make_label_pair(const Label& one, const Label& two);

class classifier {
public:
    classifier() {}
    virtual ~classifier() {}

    virtual Label classify(const FeatureVector& observation) = 0;
    static Label empty;
};

// The SVM class implements the Support Vector Machine
// discriminant function.  Instances are constructed with
// a vector of class labels (+1 or -1), a vector of dual
// coefficients, a vector of observations, and a bias value.
//
// The class SmoTrainer is responsible for determining the dual
// coefficients and bias value.
//
class SVM : public classifier {
public:
    SVM(const std::vector<double>& yy, const std::vector<double>& aa, const LabeledObservationVector& oo, double bb, const NumericClassToLabel& mm) :
        y(yy), a(aa), x(oo), b(bb), discriminantToLabel(mm) {}
    ~SVM() {}

    // the classify method should accept a list of observations?
    int discriminant(const FeatureVector& observation);
    Label classify(const FeatureVector& observation);
    double score(const LabeledObservationVector&);

private:
    const std::vector<double> y;
    const std::vector<double> a;
    const LabeledObservationVector& x;
    const double b;
    NumericClassToLabel discriminantToLabel; // trouble if this is declared const....
};


class MultiClassSVM : public classifier {
public:
	MultiClassSVM(std::vector<SVM*>);
	~MultiClassSVM();

    // the classify method should accept a list of observations
    Label classify(const FeatureVector& observation);
    double score(const LabeledObservationVector&);
private:
    const std::vector<SVM*> twoClassSvmList;
};


// the SmoTrainer trains a 2-class SVM
class SmoTrainer {
public:
	SmoTrainer();
	~SmoTrainer();

	double getC()       { return C; }
    void setC(double C) { this->C = C; }

    SVM* train(const LabeledObservationVector&);
    void assignNumericLabels(std::vector<double>&, const LabeledObservationVector&, NumericClassToLabel&);
    void elementwise_multiply(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c) {
        std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::multiplies<double>());
    }

private:
    double C = 1.0;
};


class OneVsOneMultiClassSvmTrainer {
public:
    //OneVsOneMultiClassSvmTrainer(const ObservationVector& observations, const LabelVector& observationLabels);
    OneVsOneMultiClassSvmTrainer(const LabeledObservationVector& observations);
    ~OneVsOneMultiClassSvmTrainer();

    // training requires splitting the data in to training and testing sets
    // for now, this class will take responsibility for splitting
    // return a pointer to MultiClassSVM
    MultiClassSVM* train();
    const LabelSet& getLabelSet() { return labelSet; }
    const LabeledObservationVector& getLabeledObservations() { return labeledObservations; }
    const LabelPairSet& getLabelPairSet() { return labelPairSet; }

    static void buildLabelSet(LabelSet&, const LabeledObservationVector&);
    static void buildLabelToLabeledObservationVector(LabelToLabeledObservationVector&, const LabeledObservationVector&);
    static void buildLabelPairSet(LabelPairSet&, const LabeledObservationVector&);
    static void appendTrainingAndTestingData(Label, const LabeledObservationVector&, LabeledObservationVector&, LabeledObservationVector&);
    static void standardizeObservations(const LabeledObservationVector&);


private:
    const LabeledObservationVector& labeledObservations;
    //const ObservationVector& observations;
    //const LabelVector& observationLabels;
    LabelSet labelSet;
    LabelToLabeledObservationVector labelToLabeledObservationVector;
    LabelPairSet labelPairSet;
};

// KFoldLabeledObservationDivider is used in cross validation to generate
// training and testing data sets of labeled observations.  The labels will
// be distributed in equal proportions, as much as possible.
//
// Consider a data set with 100 observations from five classes.  Also, let
// each class have 20 observations.  If we want to do 10-fold cross validation
// then training sets should have 90 observations and test sets should have
// 10 observations.  A training set should have approximately equal representation
// from each class, as should the test sets.
//
// An instance of KFoldLabeledObservationDivider will generate training and test
// sets within a for loop like this:
//
//   KFoldLabeledObservationDivider X(10, allLabeledObservations);
//   for (X.start(); !X.end(); X.next()) {
//        const LabeledObservationVector& trainingData = X.getTrainingData();
//        const LabeledObservationVector& testingData  = X.getTestingData();
//        // do cross validation on one fold
//   }
class KFoldLabeledObservationsDivider {
public:
    KFoldLabeledObservationsDivider(int _K, const LabeledObservationVector& l) : K(_K) {
        OneVsOneMultiClassSvmTrainer::buildLabelToLabeledObservationVector(labelToLabeledObservationVector, l);
    }
    ~KFoldLabeledObservationsDivider() {}

    // argument labelContainer holds the labels we will work with
    // would it be better to restrict the KFoldLabeledObservationsDivider at construction?
    void start(const LabelVector& labelContainer) {
        labelVector.clear();
        labelVector.assign(labelContainer.begin(), labelContainer.end());
        k = 0;
        next();
    }

    bool end() {
        return k > K;
    }

    void next() {
        trainingData.clear();
        testingData.clear();
        for (LabelVector::const_iterator label = labelVector.begin(); label != labelVector.end(); label++ ) {
            appendKthFold(k, K, labelToLabeledObservationVector[*label], trainingData, testingData);
        }
        k++;
    }

    int getFoldNumber() { return k; }
    const LabeledObservationVector& getTrainingData() { return trainingData; }
    const LabeledObservationVector& getTestingData() { return testingData; }

    static void appendKthFold(int k, int K, const LabeledObservationVector& x, LabeledObservationVector& trainingData, LabeledObservationVector& testingData) {
        for ( int i = 0; i < x.size(); i++) {
            if ( (i % K) == k) {
                testingData.push_back(x[i]);
            }
            else {
                trainingData.push_back(x[i]);
            }
        }
    }

private:
    int K;
    int k;
    LabelVector labelVector;
    LabelToLabeledObservationVector labelToLabeledObservationVector;
    LabeledObservationVector trainingData;
    LabeledObservationVector testingData;
};

#endif /* svm_hpp_ */
