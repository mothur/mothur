//
//  communitytype.h
//  Mothur
//
//  Created by SarahsWork on 12/3/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#ifndef Mothur_communitytype_h
#define Mothur_communitytype_h

#define EPSILON numeric_limits<double>::epsilon()


#include "mothurout.h"
#include "linearalgebra.h"
/**************************************************************************************************/

class CommunityTypeFinder {
    
public:
	CommunityTypeFinder(){	m = MothurOut::getInstance();  }
	virtual ~CommunityTypeFinder(){};
    
    virtual void printZMatrix(string, vector<string>);
    virtual void printRelAbund(string, vector<string>);
    virtual void printFitData(ofstream&) {}
    virtual void printFitData(ostream&, double) {}
    virtual void printSilData(ofstream&, double, vector<double>);
    virtual void printSilData(ostream&, double, vector<double>);
    
    virtual double getNLL()     {    return currNLL;        }
    virtual double getAIC()     {    return aic;            }
    virtual double getBIC()     {    return bic;            }
    virtual double getLogDet()  {    return logDeterminant; }
    virtual double getLaplace() {    return laplace;        }
    
    virtual double calcCHIndex(vector< vector< double> >) {return 0;}  //Calinski-Harabasz
    virtual vector<double> calcSilhouettes(vector< vector< double> >) {  vector<double> s; return s; } //if none provided by child class


protected:
    
    int findkMeans();
    vector<vector<double> > getHessian();
    double psi1(double);
    double psi(double);
    double cheb_eval(const double[], int, double);

    
	MothurOut* m;
    vector<vector<double> > zMatrix;
    vector<vector<double> > lambdaMatrix;
    vector<vector<double> > error;
    vector<vector<int> > countMatrix;
    vector<double> weights;



    int numPartitions;
    int numSamples;
    int numOTUs;
    int currentPartition;
    
    double currNLL, aic, bic, logDeterminant, laplace;
     
	
};

/**************************************************************************************************/




#endif
