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
	int tred2(vector<vector<double> >&, vector<double>&, vector<double>&);
	int qtli(vector<double>&, vector<double>&, vector<vector<double> >&);
	vector< vector<double> > calculateEuclidianDistance(vector<vector<double> >&, int);
	double calcPearson(vector<vector<double> >&, vector<vector<double> >&);
	
private:
	MothurOut* m;
	
	double pythag(double, double);
};

#endif

