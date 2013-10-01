#ifndef LINEARALGEBRA
#define LINEARALGEBRA

/*
 *  linearalgebra.h
 *  mothur
 *
 *  Created by westcott on 1/7/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "mothurout.h"


class LinearAlgebra {
	
public:
	LinearAlgebra() { m = MothurOut::getInstance(); }
	~LinearAlgebra() {}
	
	vector<vector<double> > matrix_mult(vector<vector<double> >, vector<vector<double> >);
    vector<vector<double> >transpose(vector<vector<double> >);
	void recenter(double, vector<vector<double> >, vector<vector<double> >&);
	//eigenvectors
    int tred2(vector<vector<double> >&, vector<double>&, vector<double>&);
	int qtli(vector<double>&, vector<double>&, vector<vector<double> >&);
    
	vector< vector<double> > calculateEuclidianDistance(vector<vector<double> >&, int); //pass in axes and number of dimensions
	vector< vector<double> > calculateEuclidianDistance(vector<vector<double> >&); //pass in axes
	vector<vector<double> > getObservedEuclideanDistance(vector<vector<double> >&);
	double calcPearson(vector<vector<double> >&, vector<vector<double> >&);
	double calcSpearman(vector<vector<double> >&, vector<vector<double> >&);
	double calcKendall(vector<vector<double> >&, vector<vector<double> >&);
    double calcKruskalWallis(vector<spearmanRank>&, double&);
    double calcWilcoxon(vector<double>&, vector<double>&, double&);
	
	double calcPearson(vector<double>&, vector<double>&, double&);
	double calcSpearman(vector<double>&, vector<double>&, double&);
	double calcKendall(vector<double>&, vector<double>&, double&);
    
	double calcSpearmanSig(double, double, double, double); //length, f^3 - f where f is the number of ties in x, f^3 - f where f is the number of ties in y, sum of squared diffs in ranks. - designed to find the sif of one score.
    double calcPearsonSig(double, double); //length, coeff.
    double calcKendallSig(double, double); //length, coeff.
    
    vector<double> solveEquations(vector<vector<double> >, vector<double>);
    vector<float> solveEquations(vector<vector<float> >, vector<float>);
    vector<vector<double> > getInverse(vector<vector<double> >);
    double choose(double, double);
    double normalvariate(double mu, double sigma);
    vector< vector<double> > lda(vector< vector<double> >& a, vector<string> groups, vector< vector<double> >& means, bool&); //Linear discriminant analysis - a is [features][valuesFromGroups] groups indicates which group each sampling comes from. For example if groups = early, late, mid, early, early. a[0][0] = value for feature0 from groupEarly.
    int svd(vector< vector<double> >& a, vector<double>& w, vector< vector<double> >& v); //Singular value decomposition
private:
	MothurOut* m;
	
	double pythag(double, double);
    double betacf(const double, const double, const double);
    double betai(const double, const double, const double);
    double gammln(const double);
    double gammq(const double, const double);
    double gser(double&, const double, const double, double&);
    double gcf(double&, const double, const double, double&);
    double erfcc(double);
    double gammp(const double, const double);
    double pnorm(double x);
    
    double ran0(int&); //for testing 
    double ran1(int&); //for testing
    double ran2(int&); //for testing
    double ran3(int&); //for testing
    double ran4(int&); //for testing
    void psdes(unsigned long &, unsigned long &); //for testing
    
    void ludcmp(vector<vector<double> >&, vector<int>&, double&);
    void lubksb(vector<vector<double> >&, vector<int>&, vector<double>&);
    
    void ludcmp(vector<vector<float> >&, vector<int>&, float&);
    void lubksb(vector<vector<float> >&, vector<int>&, vector<float>&);
    
};

#endif

