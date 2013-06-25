//
//  svm.hpp
//  support vector machine
//
//  Created by Joshua Lynch on 6/19/2013.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#ifndef svm_hpp_
#define svm_hpp_

#include <cmath>
#include <string>
#include <map>
#include <set>
#include <vector>

typedef std::vector<double> FeatureVector;
typedef std::vector<FeatureVector*> ObservationVector;
typedef std::string Label;
typedef std::vector<std::string> LabelVector;
typedef std::set<std::string> LabelSet;
typedef std::pair<Label,Label> LabelPair;
typedef std::set<LabelPair> LabelPairSet;
typedef std::map<Label, ObservationVector> LabeledObservations;


class classifier {
public:
	classifier() {}
	virtual ~classifier() {}

	// better to return something other than int?
	// how to represent classes?
	// int will work but ultimately a human-readable class representation is neeeded
	virtual int classify(double observations) = 0;
};

class SVM : public classifier {
private:
	// need a set of weights
	// need dual coefficients??


public:
	SVM();
	~SVM();

	// stub
	// the classify method should accept a list of observations
	int classify(double observations) { return 0; }
	double score(const ObservationVector& twoClassObservationVector, const LabelVector& twoClassLabelVector) { return 0.0; }
};


class MultiClassSVM : public classifier {
private:
	// need a set of two-class SVM classifiers

public:
	MultiClassSVM();
	~MultiClassSVM();

	// stub
	// the classify method should accept a list of observations
	void classify(const std::vector<double> observation, int& classification) {}

	void classify(const std::vector<std::vector<double> > observations, std::vector<int> classifications) {}
};


// the SmoTrainer trains a 2-class SVM
class SmoTrainer {
public:
	SmoTrainer();
	~SmoTrainer();

	double getC()       { return this->C; }
    void setC(double C) { this->C = C; }

    SVM* train(const ObservationVector& twoClassObservationVector, const LabelVector& twoClassLabelVector) {return NULL;}

private:
    double C = nan("");
};


class OneVsOneMultiClassSvmTrainer {
public:
    OneVsOneMultiClassSvmTrainer(const ObservationVector& observations, const LabelVector& observationLabels);
    ~OneVsOneMultiClassSvmTrainer();

    // training requires splitting the data in to training and testing sets
    // for now, this class will take responsibility for splitting
    // return a pointer to MultiClassSVM
    SVM* train(const ObservationVector& observations, const LabelVector& observationLabels);
    const LabelSet& getLabelSet() { return labelSet; }
    void getLabeledObservations(LabeledObservations&, const ObservationVector& observations, const LabelVector& observationLabels);
    void getLabelPairSet(LabelPairSet& labelPairSet, const LabelSet& labelSet);

private:
    const ObservationVector& observations;
    const LabelVector& observationLabels;
    LabelSet labelSet;
    LabeledObservations labeledObservations;
};

#endif /* svm_hpp_ */
