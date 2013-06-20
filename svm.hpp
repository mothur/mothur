//
//  svm.hpp
//  support vector machine
//
//  Created by Joshua Lynch on 6/19/2013.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#ifndef svm_hpp_
#define svm_hpp_

#include <string>
//#include <vector>


class classifier {
public:
	classifier() {}
	virtual ~classifier();

	//virtual int classify(double observations) = 0;
};

class SVM {
private:
	// need a set of weights
	// need dual coefficients??

public:
	SVM();
	virtual ~SVM();

	// stub
	// the classify method should accept a list of observations
	//int classify(double observations);
};


class MultiClassSVM {
private:
	// need a set of SVM

public:
	MultiClassSVM();
	virtual ~MultiClassSVM();

	// stub
	// the classify method should accept a list of observations
	//int classify(double observations);
};


// the SmoTrainer trains a 2-class SVM
class SmoTrainer {
public:
	SmoTrainer();
	virtual ~SmoTrainer();
};


class OneVsOneMultiClassSvmTrainer {
public:
	OneVsOneMultiClassSvmTrainer();
	virtual ~OneVsOneMultiClassSvmTrainer();

    // need to specify at least observations and labels
	// optionally specify cross-validation details
	//SVM* train(const std::vector < std::vector<double> > dataSet, const std::vector<std::string&> labels);
};

#endif /* svm_hpp_ */
