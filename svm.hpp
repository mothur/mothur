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
#include <deque>
#include <exception>
#include <list>
#include <map>
#include <set>
#include <stack>
#include <string>
#include <sstream>
#include "mothurout.h"


// For the purpose of training a support vector machine
// we need to calculate a dot product between two feature
// vectors.  In general these feature vectors are not
// restricted to lists of doubles, but in this implementation
// feature vectors (or 'observations' as they will be called from here on)
// will be vectors of doubles.
typedef vector<double> Observation;

/*
class Observation {
public:
    Observation() {}
    ~Observation() {}

private:
    vector<double> obs;
};
*/

// A dataset is a collection of labeled observations.
// The ObservationVector typedef is a vector
// of pointers to ObservationVectors.  Pointers are used here since
// datasets will be rearranged many times during cross validation.
// Using pointers to Observations makes copying the elements of
// an ObservationVector cheap.
typedef vector<Observation*> ObservationVector;

// Training a support vector machine requires labeled data.  The
// Label typedef defines what will constitute a class 'label' in
// this implementation.
typedef string Label;
typedef vector<Label> LabelVector;
typedef set<Label> LabelSet;

// Pairs of class labels are important because a support vector machine
// can only learn two classes of data.  The LabelPair typedef is a vector
// even though a pair might seem more natural, but it is useful to
// iterate over the pair.
typedef vector<Label> LabelPair;
LabelPair buildLabelPair(const Label& one, const Label& two);

// Learning to classify a dataset with more than two classes requires
// training a separate support vector machine for each pair of classes.
// The LabelPairSet typedef defines a container for the collection of
// all unique label pairs for a set of data.
typedef set<LabelPair> LabelPairSet;

// A dataset is a set of observations with associated labels.  The
// LabeledObservation typedef is a label-observation pair intended to
// hold one observation and its corresponding label.  Using a pointer
// to Observation makes these objects cheap to copy.
//typedef pair<Label, Observation*> LabeledObservation;

// This is a refactoring of the original LabeledObservation typedef.
// The original typedef has been promoted to a class in order to add
// at least one additional member variable, int datasetIndex, which
// will be used to implement kernel function optimizations.
class LabeledObservation {
public:
    LabeledObservation(int _datasetIndex, Label _label, Observation* _o) : datasetIndex(_datasetIndex), first(_label), second(_o) {}
    ~LabeledObservation() {}

    void removeFeatureAtIndex(int n) {
        int m = 0;
        Observation::iterator i = second->begin();
        while ( m < n ) {
            i++;
            m++;
        }
        second->erase(i);
    }

    int getDatasetIndex()         const   { return datasetIndex; }
    Label getLabel()              const   { return first; }
    Observation* getObservation() const   { return second; }

//private:
    int datasetIndex;
    Label first;
    Observation* second;
};



// A LabeledObservationVector is a container for an entire dataset (or a
// subset of an entire dataset).
typedef vector<LabeledObservation> LabeledObservationVector;
void buildLabelSet(LabelSet&, const LabeledObservationVector&);


double getMinimumFeatureValueForObservation(Observation::size_type featureIndex, LabeledObservationVector& observations);
double getMaximumFeatureValueForObservation(Observation::size_type featureIndex, LabeledObservationVector& observations);


void transformZeroOne(LabeledObservationVector&);
void transformZeroMeanUnitVariance(LabeledObservationVector&);


class Feature {
public:
    Feature(int i, const string& l) : index(i), label(l) {}
    Feature(const Feature& f) : index(f.index), label(f.label) {}
    ~Feature() {}

    int getFeatureIndex() const { return index; }
    void setFeatureIndex(int i) { index = i; }
    string getFeatureLabel() const { return label; }

private:
    int index;
    string label;
};

typedef list<Feature> FeatureList;
typedef vector<Feature> FeatureVector;

// might make sense for this to be a member function of SvmDataset
FeatureVector applyStdThreshold(double, LabeledObservationVector&, FeatureVector&);


//  A RankedFeature is just a Feature and a its associated 'rank', where
//  rank is the SVM-RFE iteration during which the feature was eliminated.
//  If the SVM-RFE method eliminates multiple features in an iteration
//  then some features will have the same rank.
class RankedFeature {
public:
    RankedFeature(const Feature& f, int r) : feature(f), rank(r) {}
    ~RankedFeature() {}

    Feature getFeature() const { return feature; }
    int getRank() const { return rank; }

private:
    Feature feature;
    int rank;
};

typedef list<RankedFeature> RankedFeatureList;


// The SvmDataset class encapsulates labeled observations and feature information.
// All data required to train SVMs is found in SvmDataset.
class SvmDataset {
public:
    SvmDataset(const LabeledObservationVector& v, const FeatureVector& f) : labeledObservationVector(v), featureVector(f) {}
    ~SvmDataset() {}

    LabeledObservationVector& getLabeledObservationVector() { return labeledObservationVector; }
    FeatureVector& getFeatureVector() { return featureVector; }

    void removeFeature(const Feature feature) {

    }

private:
    LabeledObservationVector labeledObservationVector;
    FeatureVector featureVector;
};


//
//  0 - print no optional information (quiet)
//  1 - print minimum optional information (info)
//  2 - print a little more optional information (debug)
//  3 - print the maximum amount of optional information (trace)
//
class OutputFilter {
public:
    OutputFilter(int v) : verbosity(v) {}
    OutputFilter(const OutputFilter& of) : verbosity(of.verbosity) {}
    ~OutputFilter() {}

    bool info()  const { return verbosity >= INFO; }
    bool debug() const { return verbosity >= DEBUG; }
    bool trace() const { return verbosity >= TRACE; }

    static const int QUIET;
    static const int INFO;
    static const int DEBUG;
    static const int TRACE;

private:
    const int verbosity;
};


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
typedef map<Label, LabeledObservationVector> LabelToLabeledObservationVector;
void buildLabelToLabeledObservationVector(LabelToLabeledObservationVector&, const LabeledObservationVector&);

// A support vector machine uses +1 and -1 in calculations to represent
// the two classes of data it is trained to distinguish.  The NumericClassToLabel
// container is used to record the labels associated with these integers.
// For a dataset with labels "blue" and "green" a NumericClassToLabel map looks like
//     1, "blue"
//    -1, "green"
typedef map<int, Label> NumericClassToLabel;
void buildNumericClassToLabelMap(LabelPair);

typedef double Parameter;
typedef string ParameterName;
typedef vector<double> ParameterRange;
typedef map<ParameterName, ParameterRange> ParameterRangeMap;

typedef map<string, ParameterRangeMap> KernelParameterRangeMap;
void getDefaultKernelParameterRangeMap(KernelParameterRangeMap& kernelParameterRangeMap);

typedef map<ParameterName, Parameter> ParameterMap;
typedef vector<ParameterMap> ParameterMapVector;
typedef stack<Parameter> ParameterStack;

class ParameterSetBuilder {
public:
    // If the argument ParameterRangeMap looks like this:
    //     { "a" : [1.0, 2.0], "b" : [-1.0, 1.0], "c" : [0.5, 0.6] }
    // then the list of parameter sets looks like this:
    //     [ {"a":1.0, "b":-1.0, "c":0.5},
    //       {"a":1.0, "b":-1.0, "c":0.6},
    //       {"a":1.0, "b": 1.0, "c":0.5},
    //       {"a":1.0, "b": 1.0, "c":0.6},
    //       {"a":2.0, "b":-1.0, "c":0.5},
    //       {"a":2.0, "b":-1.0, "c":0.6},
    //       {"a":2.0, "b": 1.0, "c":0.5},
    //       {"a":2.0, "b": 1.0, "c":0.6},
    //     ]
    ParameterSetBuilder(const ParameterRangeMap& parameterRangeMap) {
        // a small step toward quieting down this code
        bool verbose = false;

        stack<pair<ParameterName, ParameterStack> > stackOfParameterRanges;
        stack<pair<ParameterName, ParameterStack> > stackOfEmptyParameterRanges;
        ParameterMap nextParameterSet;
        int parameterSetCount = 1;
        for ( ParameterRangeMap::const_iterator i = parameterRangeMap.begin(); i != parameterRangeMap.end(); i++ ) {
            parameterSetCount *= i->second.size();
            ParameterName parameterName = i->first;
            ParameterStack emptyParameterStack;
            stackOfEmptyParameterRanges.push(make_pair(parameterName, emptyParameterStack));
        }
        // get started
        for ( int n = 0; n < parameterSetCount; n++ ) {

            if (verbose) m->mothurOut("n = " + toString(n) ); m->mothurOutEndLine();

            // pull empty stacks off until there are no empty stacks
            while ( stackOfParameterRanges.size() > 0 and stackOfParameterRanges.top().second.size() == 0 ) {

                if (verbose) m->mothurOut("  empty parameter range: " + stackOfParameterRanges.top().first); m->mothurOutEndLine();

                stackOfEmptyParameterRanges.push(stackOfParameterRanges.top());
                stackOfParameterRanges.pop();
            }

            // move to the next value for the parameter at the top of the stackOfParameterRanges
            if ( stackOfParameterRanges.size() > 0 ) {
                if (verbose) {
                    m->mothurOut( "  moving to next value for parameter " + toString(stackOfParameterRanges.top().first) ); m->mothurOutEndLine();
                    m->mothurOut( "    next value is "  + toString(stackOfParameterRanges.top().second.top()) ); m->mothurOutEndLine();
                }
                ParameterName parameterName = stackOfParameterRanges.top().first;
                nextParameterSet[parameterName] = stackOfParameterRanges.top().second.top();
                stackOfParameterRanges.top().second.pop();
            }
            if (verbose) m->mothurOut( "stack of empty parameter ranges has size " + toString(stackOfEmptyParameterRanges.size() ) ); m->mothurOutEndLine();
            // reset each parameter range that has been exhausted
            while ( stackOfEmptyParameterRanges.size() > 0 ) {
                ParameterName parameterName = stackOfEmptyParameterRanges.top().first;
                if (verbose) m->mothurOut( "  reseting range for parameter " + toString(stackOfEmptyParameterRanges.top().first) ); m->mothurOutEndLine();
                stackOfParameterRanges.push(stackOfEmptyParameterRanges.top());
                stackOfEmptyParameterRanges.pop();
                const ParameterRange& parameterRange = parameterRangeMap.find(parameterName)->second;
                // it is nice to have the parameters used in order smallest to largest
                // so that we choose the smallest in ties
                // but we will not enforce this so users can specify parameters in the order they like
                // this loop will use parameters in the order they are found in the parameter range
                for (ParameterRange::const_reverse_iterator i = parameterRange.rbegin(); i != parameterRange.rend(); i++ ) {
                    stackOfParameterRanges.top().second.push(*i);
                }
                nextParameterSet[parameterName] = stackOfParameterRanges.top().second.top();
                stackOfParameterRanges.top().second.pop();
            }
            parameterSetVector.push_back(nextParameterSet);
            // print out the next parameter set
            if (verbose) {
                for (ParameterMap::iterator p = nextParameterSet.begin(); p != nextParameterSet.end(); p++) {
                    m->mothurOut(toString(p->first) + " : " + toString(p->second) ); m->mothurOutEndLine();
                }
            }
        }
    }
    ~ParameterSetBuilder() {}

    const ParameterMapVector& getParameterSetList() { return parameterSetVector; }

private:
    ParameterMapVector parameterSetVector;
    MothurOut* m;
};


class RowCache {
public:
    RowCache(int d) : cache(d, NULL) {}
    virtual ~RowCache() {
        for (int i = 0; i < cache.size(); i++) {
            if ( !rowNotCached(i) ) {
                delete cache[i];
            }
        }
    }

    double getCachedValue(int i, int j) {
        if ( rowNotCached(i) ) {
            createRow(i);
        }
        return cache.at(i)->at(j);
    }

    void createRow(int i) {
        cache[i] = new vector<double>(cache.size(), numeric_limits<double>::signaling_NaN());
        for ( int v = 0; v < cache.size(); v++ ) {
            cache.at(i)->at(v) = calculateValueForCache(i, v);
        }
    }

    bool rowNotCached(int i) {
        return cache[i] == NULL;
    }

    virtual double calculateValueForCache(int, int) = 0;

private:
    vector<vector<double>* > cache;
};


class InnerProductRowCache : public RowCache {
public:
    InnerProductRowCache(const LabeledObservationVector& _obs) : obs(_obs), RowCache(_obs.size()) {}
    virtual ~InnerProductRowCache() {}

    double getInnerProduct(const LabeledObservation& obs_i, const LabeledObservation& obs_j) {
        return getCachedValue(
            obs_i.datasetIndex,
            obs_j.datasetIndex
        );
    }

    double calculateValueForCache(int i, int j) {
        return inner_product(obs[i].second->begin(), obs[i].second->end(), obs[j].second->begin(), 0.0);
    }

private:
    const LabeledObservationVector& obs;
};


// The KernelFunction class caches a partial kernel value that does not depend on kernel parameters.
class KernelFunction {
public:
    //KernelFunction(const LabeledObservationVector& _obs, InnerProductCache& _ipc) : obs(_obs), innerProductRowCache(_ipc) {}

    KernelFunction(const LabeledObservationVector& _obs) :
        obs(_obs),
        cache(_obs.size(), NULL) {}

    virtual ~KernelFunction() {
        for (int i = 0; i < cache.size(); i++) {
            if ( !rowNotCached(i) ) {
                delete cache[i];
            }
        }
    }

    virtual double similarity(const LabeledObservation&, const LabeledObservation&) = 0;
    virtual void setParameters(const ParameterMap&) = 0;
    virtual void getDefaultParameterRanges(ParameterRangeMap&) = 0;

    virtual double calculateParameterFreeSimilarity(const LabeledObservation&, const LabeledObservation&) = 0;

    double getCachedParameterFreeSimilarity(const LabeledObservation& obs_i, const LabeledObservation& obs_j) {
        const int i = obs_i.datasetIndex;
        const int j = obs_j.datasetIndex;

        if ( rowNotCached(i) ) {
            cache[i] = new vector<double>(obs.size(), numeric_limits<double>::signaling_NaN());
            for ( int v = 0; v < obs.size(); v++ ) {
                cache.at(i)->at(v) = calculateParameterFreeSimilarity(obs[i], obs[v]);
            }
        }
        return cache.at(i)->at(j);
    }

    bool rowNotCached(int i) {
        return cache[i] == NULL;
    }

private:
    const LabeledObservationVector& obs;
    //vector<vector<double> > cache;
    vector<vector<double>* > cache;
    //InnerProductRowCache& innerProductCache;
};


class LinearKernelFunction : public KernelFunction {
public:
    // parameters must be set before using a KernelFunction is used
    LinearKernelFunction(const LabeledObservationVector& _obs) : KernelFunction(_obs), constant(0.0) {}
    ~LinearKernelFunction() {}

    double similarity(const LabeledObservation& i, const LabeledObservation& j) {
        return getCachedParameterFreeSimilarity(i, j) + constant;
    }

    double calculateParameterFreeSimilarity(const LabeledObservation& i, const LabeledObservation& j) {
        return inner_product(i.second->begin(), i.second->end(), j.second->begin(), 0.0);
    }

    double getConstant() { return constant; }
    void setConstant(double c) { constant = c; }

    void setParameters(const ParameterMap& p) {
        setConstant(p.find(MapKey_Constant)->second);
    };

    void getDefaultParameterRanges(ParameterRangeMap& p) {
        p[MapKey_Constant] = defaultConstantRange;
    }

    static const string MapKey;
    static const string MapKey_Constant;
    static const ParameterRange defaultConstantRange;

private:
    double constant;
};


class RbfKernelFunction : public KernelFunction {
public:
    // parameters must be set before a KernelFunction is used
    RbfKernelFunction(const LabeledObservationVector& _obs) : KernelFunction(_obs), gamma(0.0) {}
    ~RbfKernelFunction() {}

    double similarity(const LabeledObservation& i, const LabeledObservation& j) {
        //double sumOfSquaredDifs = 0.0;
        //for (int n = 0; n < i.second->size(); n++) {
        //    sumOfSquaredDifs += pow((i.second->at(n) - j.second->at(n)), 2.0);
        //}
        return gamma * getCachedParameterFreeSimilarity(i, j);
    }

    double calculateParameterFreeSimilarity(const LabeledObservation& i, const LabeledObservation& j) {
        //double sumOfSquaredDifs = 0.0;
        //for (int n = 0; n < i.second->size(); n++) {
        //    sumOfSquaredDifs += pow((i.second->at(n) - j.second->at(n)), 2.0);
        //}
        double sumOfSquaredDifs =
                      inner_product(i.second->begin(), i.second->end(), i.second->begin(), 0.0)
              - 2.0 * inner_product(i.second->begin(), i.second->end(), j.second->begin(), 0.0)
              +       inner_product(j.second->begin(), j.second->end(), j.second->begin(), 0.0);
        return exp(sqrt(sumOfSquaredDifs));
    }

    double getGamma()       { return gamma; }
    void setGamma(double g) { gamma = g; }

    void setParameters(const ParameterMap& p) {
        setGamma(p.find(MapKey_Gamma)->second);
    }

    void getDefaultParameterRanges(ParameterRangeMap& p) {
        p[MapKey_Gamma] = defaultGammaRange;
    }

    static const string MapKey;
    static const string MapKey_Gamma;

    static const ParameterRange defaultGammaRange;

private:
    double gamma;
};


class PolynomialKernelFunction : public KernelFunction {
public:
    // parameters must be set before using a KernelFunction is used
    PolynomialKernelFunction(const LabeledObservationVector& _obs) : KernelFunction(_obs), c(0.0), gamma(0.0), d(0) {}
    ~PolynomialKernelFunction() {}

    double similarity(const LabeledObservation& i, const LabeledObservation& j) {
        return pow((gamma * getCachedParameterFreeSimilarity(i, j) + c), d);
        //return pow(inner_product(i.second->begin(), i.second->end(), j.second->begin(), c), d);
    }

    double calculateParameterFreeSimilarity(const LabeledObservation& i, const LabeledObservation& j) {
        return inner_product(i.second->begin(), i.second->end(), j.second->begin(), 0.0);
    }

    void setParameters(const ParameterMap& p) {
        c = p.find(MapKey_Constant)->second;
        gamma = p.find(MapKey_Coefficient)->second;
        d = int(p.find(MapKey_Degree)->second);
    }

    void getDefaultParameterRanges(ParameterRangeMap& p) {
        p[MapKey_Constant] = defaultConstantRange;
        p[MapKey_Coefficient] = defaultCoefficientRange;
        p[MapKey_Degree] = defaultDegreeRange;
    }

    static const string MapKey;
    static const string MapKey_Constant;
    static const string MapKey_Coefficient;
    static const string MapKey_Degree;

    static const ParameterRange defaultConstantRange;
    static const ParameterRange defaultCoefficientRange;
    static const ParameterRange defaultDegreeRange;

private:
    double c;
    double gamma;
    int d;
};


class SigmoidKernelFunction : public KernelFunction {
public:
    // parameters must be set before using a KernelFunction is used
    SigmoidKernelFunction(const LabeledObservationVector& _obs) : KernelFunction(_obs), alpha(0.0), c(0.0) {}
    ~SigmoidKernelFunction() {}

    double similarity(const LabeledObservation& i, const LabeledObservation& j) {
        return tanh(alpha * getCachedParameterFreeSimilarity(i, j) + c);
        //return tanh(alpha * inner_product(i.second->begin(), i.second->end(), j.second->begin(), c));
    }

    double calculateParameterFreeSimilarity(const LabeledObservation& i, const LabeledObservation& j) {
        return inner_product(i.second->begin(), i.second->end(), j.second->begin(), 0.0);
    }

    void setParameters(const ParameterMap& p) {
        alpha = p.find(MapKey_Alpha)->second;
        c = p.find(MapKey_Constant)->second;
    }

    void getDefaultParameterRanges(ParameterRangeMap& p) {
        p[MapKey_Alpha] = defaultAlphaRange;
        p[MapKey_Constant] = defaultConstantRange;
    }

    static const string MapKey;
    static const string MapKey_Alpha;
    static const string MapKey_Constant;

    static const ParameterRange defaultAlphaRange;
    static const ParameterRange defaultConstantRange;
private:
    double alpha;
    double c;
};


class KernelFactory {
public:
    static KernelFunction* getKernelFunctionForKey(string kernelFunctionKey, const LabeledObservationVector& obs) {
        if ( kernelFunctionKey == LinearKernelFunction::MapKey ) {
            return new LinearKernelFunction(obs);
        }
        else if ( kernelFunctionKey == RbfKernelFunction::MapKey ) {
            return new RbfKernelFunction(obs);
        }
        else if ( kernelFunctionKey == PolynomialKernelFunction::MapKey ) {
            return new PolynomialKernelFunction(obs);
        }
        else if ( kernelFunctionKey == SigmoidKernelFunction::MapKey ) {
            return new SigmoidKernelFunction(obs);
        }
        else {
            throw new exception();
        }
    }
};


typedef map<string, KernelFunction*> KernelFunctionMap;

// An instance of KernelFunctionFactory dynamically allocates kernel function
// instances and maintains a table of pointers to them.  This allows kernel
// function instances to be reused which improves performance since the
// kernel values do not have to be recalculated as often.  A KernelFunctionFactory
// maintains an inner product cache used by the KernelFunctions it builds.
class KernelFunctionFactory {
public:
    KernelFunctionFactory(const LabeledObservationVector& _obs) : obs(_obs) {}
    ~KernelFunctionFactory() {
        for ( KernelFunctionMap::iterator i = kernelFunctionTable.begin(); i != kernelFunctionTable.end(); i++ ) {
            delete i->second;
        }
    }

    KernelFunction& getKernelFunctionForKey(string kernelFunctionKey) {
        if ( kernelFunctionTable.count(kernelFunctionKey) == 0 ) {
            kernelFunctionTable.insert(
                make_pair(
                    kernelFunctionKey,
                    KernelFactory::getKernelFunctionForKey(kernelFunctionKey, obs)
                )
            );
        }
        return *kernelFunctionTable[kernelFunctionKey];
    }

private:
    const LabeledObservationVector& obs;
    KernelFunctionMap kernelFunctionTable;
    //InnerProductCache innerProductCache;
};


class KernelFunctionCache {
public:
    KernelFunctionCache(KernelFunction& _k, const LabeledObservationVector& _obs) :
        k(_k), obs(_obs),
        cache(_obs.size(), NULL) {}
    ~KernelFunctionCache() {
        //cout << "deleting KernelFunctionCache cache" << endl;
        for (int i = 0; i < cache.size(); i++) {
            if ( !rowNotCached(i) ) {
                delete cache[i];
            }
        }
    }

    double similarity(const LabeledObservation& obs_i, const LabeledObservation& obs_j) {
        const int i = obs_i.datasetIndex;
        const int j = obs_j.datasetIndex;
        // if the first element of row i is NaN then calculate all elements for row i
        if ( rowNotCached(i) ) {
            cache[i] = new vector<double>(obs.size(), numeric_limits<double>::signaling_NaN());
            for ( int v = 0; v < obs.size(); v++ ) {
                cache.at(i)->at(v) = k.similarity(
                    obs[i],
                    obs[v]
                );
            }
        }
        return cache.at(i)->at(j);
    }

    bool rowNotCached(int i) {
        return cache[i] == NULL;
    }

private:
    KernelFunction& k;
    const LabeledObservationVector& obs;
    //vector<vector<double> > cache;
    vector<vector<double>* > cache;
};


// The SVM class implements the Support Vector Machine
// discriminant function.  Instances are constructed with
// a vector of class labels (+1.0 or -1.0), a vector of dual
// coefficients, a vector of observations, and a bias value.
//
// The class SmoTrainer is responsible for determining the dual
// coefficients and bias value.
//
class SVM {
public:
    SVM(const vector<double>& yy, const vector<double>& aa, const LabeledObservationVector& oo, double bb, const NumericClassToLabel& mm) :
        y(yy), a(aa), x(oo), b(bb), discriminantToLabel(mm) {}
    ~SVM() {}

    // the classify method should accept a list of observations?
    int discriminant(const Observation&) const;
    Label classify(const Observation& observation) const {
        //return discriminantToLabel[discriminant(observation)];
        return discriminantToLabel.find(discriminant(observation))->second;
    }
    LabelVector classify(const LabeledObservationVector&) const;
    double score(const LabeledObservationVector&) const;

    NumericClassToLabel getDiscriminantToLabel() const { return discriminantToLabel; }
    LabelPair getLabelPair() const { return buildLabelPair(discriminantToLabel.find(1)->second, discriminantToLabel.find(-1)->second); }

public:
    // y holds the numeric class: +1.0 or -1.0
    const vector<double> y;
    // a holds the optimal dual coefficients
    const vector<double> a;
    // x holds the support vectors
    const LabeledObservationVector x;
    const double b;
    const NumericClassToLabel discriminantToLabel;
};


class SvmPerformanceSummary {
public:
    SvmPerformanceSummary() {}
    // this constructor should be used by clients other than tests
    SvmPerformanceSummary(const SVM& svm, const LabeledObservationVector& actual) {
        init(svm, actual, svm.classify(actual));
    }
    // this constructor is intended for unit testing
    SvmPerformanceSummary(const SVM& svm, const LabeledObservationVector& actual, const LabelVector& predictions) {
        init(svm, actual, predictions);
    }

    Label getPositiveClassLabel() const { return positiveClassLabel; }
    Label getNegativeClassLabel() const { return negativeClassLabel; }

    double getPrecision() const { return precision; }
    double getRecall()    const { return recall; }
    double getF()         const { return f; }
    double getAccuracy()  const { return accuracy; }

private:
    void init(const SVM&, const LabeledObservationVector&, const LabelVector&);

    //const SVM& svm;

    Label positiveClassLabel;
    Label negativeClassLabel;

    double precision;
    double recall;
    double f;
    double accuracy;
};


class MultiClassSvmClassificationTie : public exception {
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

typedef vector<SVM*> SvmVector;
typedef map<LabelPair, SvmPerformanceSummary> SvmToSvmPerformanceSummary;

// Using SVM with more than two classes requires training multiple SVMs.
// The MultiClassSVM uses a vector of trained SVMs to do classification
// on data having more than two classes.
class MultiClassSVM {
public:
    MultiClassSVM(const vector<SVM*>, const LabelSet&, const SvmToSvmPerformanceSummary&, OutputFilter);
    ~MultiClassSVM();

    // the classify method should accept a list of observations
    Label classify(const Observation& observation);
    double score(const LabeledObservationVector&);

    // no need to delete these pointers
    const SvmVector& getSvmList() { return twoClassSvmList; }

    const LabelSet& getLabels() { return labelSet; }

    const SvmPerformanceSummary& getSvmPerformanceSummary(const SVM& svm) { return svmToSvmPerformanceSummary.at(svm.getLabelPair()); }

    double getAccuracy() { return accuracy; }
    void setAccuracy(const LabeledObservationVector& obs) { accuracy = score(obs); }

private:
    const SvmVector twoClassSvmList;
    const LabelSet labelSet;
    const OutputFilter outputFilter;

    double accuracy;
    MothurOut* m;

    // this is a map from label pairs to performance summaries
    SvmToSvmPerformanceSummary svmToSvmPerformanceSummary;
};


//class SvmTrainingInterruptedException : public exception {
//public:
//    SvmTrainingInterruptedException(const string& m) : message(m) {}
//    ~SvmTrainingInterruptedException() throw() {}
//    virtual const char* what() const throw() {
//        return message.c_str();
//    }

//private:
//    string message;
//};

class SmoTrainerException : public exception {
public:
    SmoTrainerException(const string& m) : message(m) {}
    ~SmoTrainerException() throw() {}
    virtual const char* what() const throw() {
        return message.c_str();
    }

private:
    string message;
};

//class ExternalSvmTrainingInterruption {
//public:
//    ExternalSvmTrainingInterruption() {}
//    virtual ~ExternalSvmTrainingInterruption() throw() {}
//    virtual bool interruptTraining() { return false; }
//};


// SmoTrainer trains a support vector machine using Sequential
// Minimal Optimization as described in the article
// "Support Vector Machine Solvers" by Bottou and Lin.
class SmoTrainer {
public:
    SmoTrainer(OutputFilter of) : outputFilter(of), C(1.0) {}

    ~SmoTrainer() {}

    double getC()       { return C; }
    void setC(double C) { this->C = C; }

    void setParameters(const ParameterMap& p) {
        C = p.find(MapKey_C)->second;
    }

    SVM* train(KernelFunctionCache&, const LabeledObservationVector&);
    void assignNumericLabels(vector<double>&, const LabeledObservationVector&, NumericClassToLabel&);
    void elementwise_multiply(vector<double>& a, vector<double>& b, vector<double>& c) {
        transform(a.begin(), a.end(), b.begin(), c.begin(), multiplies<double>());
    }

    static const string MapKey_C;
    static const ParameterRange defaultCRange;

private:
    //ExternalSvmTrainingInterruption& externalSvmTrainingInterruption;

    const OutputFilter outputFilter;

    double C;
};


// KFoldLabeledObservationDivider is used in cross validation to generate
// training and testing data sets of labeled observations.  The labels will
// be distributed in proportion to their frequency in the data, as much as possible.
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
    // initialize the k member variable to K so end() will return true if it is called before start()
    // this is not perfect protection against misuse but it's better than nothing
    KFoldLabeledObservationsDivider(int _K, const LabeledObservationVector& l) : K(_K), k(_K) {
        buildLabelToLabeledObservationVector(labelToLabeledObservationVector, l);
    }
    ~KFoldLabeledObservationsDivider() {}

    void start() {
        k = 0;
        trainingData.clear();
        testingData.clear();
        for (LabelToLabeledObservationVector::const_iterator p = labelToLabeledObservationVector.begin(); p != labelToLabeledObservationVector.end(); p++) {
            appendKthFold(k, K, p->second, trainingData, testingData);
        }
    }

    bool end() {
        return k >= K;
    }

    void next() {
        k++;
        trainingData.clear();
        testingData.clear();
        for (LabelToLabeledObservationVector::const_iterator p = labelToLabeledObservationVector.begin(); p != labelToLabeledObservationVector.end(); p++) {
            appendKthFold(k, K, p->second, trainingData, testingData);
        }
    }

    int getFoldNumber() { return k; }
    const LabeledObservationVector& getTrainingData() { return trainingData; }
    const LabeledObservationVector& getTestingData() { return testingData; }

    // Function appendKthFold takes care of partitioning the observations in x into two sets,
    // one for training and one for testing.  The argument K specifies how many folds
    // will be requested in all.  The argument k specifies which fold to return.
    // An example: let K=3, k=0, and let there be 10 observations (all having the same label)
    //     i  i%3  (i%3)==0  k=0 partition  (i%3)==1  k=1 partition  (i%3)==2  k=2 partition
    //     0   0     true      testing        false     training       false     training
    //     1   1     false     training       true      testing        false     training
    //     2   2     false     training       false     training       true      testing
    //     3   0     true      testing        false     training       false     training
    //     4   1     false     training       true      testing        false     training
    //     5   2     false     training       false     training       true      testing
    //     6   0     true      testing        false     training       false     training
    //     7   1     false     training       true      testing        false     training
    //     8   2     false     training       false     training       true      testing
    //     9   0     true      testing        false     training       false     training
    //
    static void appendKthFold(int k, int K, const LabeledObservationVector& x, LabeledObservationVector& trainingData, LabeledObservationVector& testingData) {
        //for ( int i = 0; i < x.size(); i++) {
        int i = 0;
        for (LabeledObservationVector::const_iterator xi = x.begin(); xi != x.end(); xi++) {
            if ( (i % K) == k) {
                testingData.push_back(*xi);
            }
            else {
                trainingData.push_back(*xi);
            }
            i++;
        }
    }

private:
    const int K;
    int k;
    LabelVector labelVector;
    LabelToLabeledObservationVector labelToLabeledObservationVector;
    LabeledObservationVector trainingData;
    LabeledObservationVector testingData;
};


// OneVsOneMultiClassSvmTrainer trains a support vector machine for each
// pair of labels in a set of data.
class OneVsOneMultiClassSvmTrainer {
public:
    OneVsOneMultiClassSvmTrainer(SvmDataset&, int, int, OutputFilter&);
    ~OneVsOneMultiClassSvmTrainer() {}

    MultiClassSVM* train(const KernelParameterRangeMap&);
    double trainOnKFolds(SmoTrainer&, KernelFunctionCache&, KFoldLabeledObservationsDivider&);
    const LabelSet& getLabelSet() { return labelSet; }
    const LabeledObservationVector& getLabeledObservations() { return svmDataset.getLabeledObservationVector(); }
    const LabelPairSet& getLabelPairSet() { return labelPairSet; }
    const LabeledObservationVector& getLabeledObservationVectorForLabel(const Label& label) { return labelToLabeledObservationVector[label]; }

    const OutputFilter& getOutputFilter() { return outputFilter; }

    static void buildLabelPairSet(LabelPairSet&, const LabeledObservationVector&);
    static void appendTrainingAndTestingData(Label, const LabeledObservationVector&, LabeledObservationVector&, LabeledObservationVector&);

private:

    const OutputFilter outputFilter;
    //bool verbose;

    SvmDataset& svmDataset;

    const int evaluationFoldCount;
    const int trainFoldCount;

    LabelSet labelSet;
    LabelToLabeledObservationVector labelToLabeledObservationVector;
    LabelPairSet labelPairSet;

};

// A better name for this class is MsvmRfe after MSVM-RFE described in
// "MSVM-RFE: extensions of SVM-RFE for multiclass gene selection on
// DNA microarray data", Zhou and Tuck, 2007, Bioinformatics
class SvmRfe {
public:
    SvmRfe() {}
    ~SvmRfe() {}

    RankedFeatureList getOrderedFeatureList(SvmDataset&, OneVsOneMultiClassSvmTrainer&, const ParameterRange&, const ParameterRange&);
};


#endif /* svm_hpp_ */
