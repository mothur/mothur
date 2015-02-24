#ifndef LIBSHUFF
#define LIBSHUFF

/*
 *  libshuff.h
 *  Mothur
 *
 *  Created by Pat Schloss on 4/8/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "fullmatrix.h"

class Libshuff {
	
public:
	Libshuff(FullMatrix*, int, float, float);
    virtual ~Libshuff() {}
	virtual vector<vector<double> > evaluateAll() = 0;
	virtual float evaluatePair(int,int) = 0;
	void randomizeGroups(int, int);
	void resetGroup(int);
	vector<vector<vector<double> > > getSavedMins();

protected:
	void initializeGroups(FullMatrix*);
	vector<double> getMinX(int);
	vector<double> getMinXY(int, int);
	
	vector<vector<vector<double> > > savedMins;
	
	
	FullMatrix* matrix;
	vector<int> groupSizes;
	vector<string> groupNames;
	vector<vector<int> > groups;
	vector<vector<int> > savedGroups;
	vector<double> minX;
	vector<double> minXY;
	float cutOff;
	int iters;
	float stepSize;
	
	int numGroups;
	MothurOut* m;
};

#endif
