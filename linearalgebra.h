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
	void recenter(double, vector<vector<double> >, vector<vector<double> >&);
	int tred2(vector<vector<double> >&, vector<double>&, vector<double>&);
	int qtli(vector<double>&, vector<double>&, vector<vector<double> >&);
	vector< vector<double> > calculateEuclidianDistance(vector<vector<double> >&, int); //pass in axes and number of dimensions
	vector< vector<double> > calculateEuclidianDistance(vector<vector<double> >&); //pass in axes
	double calcPearson(vector<vector<double> >&, vector<vector<double> >&);
	double calcSpearman(vector<vector<double> >&, vector<vector<double> >&);
	double calcKendall(vector<vector<double> >&, vector<vector<double> >&);
	
private:
	MothurOut* m;
	
	double pythag(double, double);
};

#endif

