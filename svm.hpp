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

typedef std::vector<double> FeatureVector;
typedef std::vector<FeatureVector*> ObservationVector;
typedef std::string Label;
typedef std::vector<std::string> LabelVector;
typedef std::set<std::string> LabelSet;
typedef std::pair<Label,Label> LabelPair;
typedef std::set<LabelPair> LabelPairSet;
typedef std::map<Label, ObservationVector*> LabeledObservations;
typedef std::map<int, Label> NumericClassToLabel;


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
    SVM(const std::vector<double>& yy, const std::vector<double>& aa, const ObservationVector& oo, double bb, const NumericClassToLabel& mm) :
        y(yy), a(aa), x(oo), b(bb), discriminantToLabel(mm) {}
    ~SVM() {}

    // the classify method should accept a list of observations?
    int discriminant(const FeatureVector& observation);
    Label classify(const FeatureVector& observation);
    double score(const ObservationVector& twoClassObservationVector, const LabelVector& twoClassLabelVector);

private:
    const std::vector<double> y;
    const std::vector<double> a;
    const ObservationVector& x;
    const double b;
    NumericClassToLabel discriminantToLabel; // trouble if this is declared const....
};


class MultiClassSVM : public classifier {
public:
	MultiClassSVM(std::vector<SVM*>);
	~MultiClassSVM();

    // the classify method should accept a list of observations
    Label classify(const FeatureVector& observation);
    double score(const ObservationVector& twoClassObservationVector, const LabelVector& twoClassLabelVector);
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

    SVM* train(const ObservationVector& twoClassObservationVector, const LabelVector& twoClassLabelVector);
    void assignNumericLabels(std::vector<double>& y, const LabelVector& labelVector, NumericClassToLabel&);
    void elementwise_multiply(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c) {
        std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::multiplies<double>());
    }

private:
    double C = 1.0;
};


class OneVsOneMultiClassSvmTrainer {
public:
    OneVsOneMultiClassSvmTrainer(const ObservationVector& observations, const LabelVector& observationLabels);
    ~OneVsOneMultiClassSvmTrainer();

    // training requires splitting the data in to training and testing sets
    // for now, this class will take responsibility for splitting
    // return a pointer to MultiClassSVM
    MultiClassSVM* train();
    const LabelSet& getLabelSet() { return labelSet; }
    LabeledObservations& getLabeledObservations() { return labeledObservations; }
    const LabelPairSet& getLabelPairSet() { return labelPairSet; }

    static void buildLabeledObservations(LabeledObservations&, const LabelSet&, const ObservationVector&, const LabelVector&);
    static void buildLabelPairSet(LabelPairSet&, const LabelSet&);
    static void appendTrainingAndTestingData(Label, const ObservationVector&, ObservationVector&, LabelVector&, ObservationVector&, LabelVector&);
    static void standardizeObservations(const ObservationVector&);


private:
    const ObservationVector& observations;
    const LabelVector& observationLabels;
    LabelSet labelSet;
    LabeledObservations labeledObservations;
    LabelPairSet labelPairSet;
};

#endif /* svm_hpp_ */
