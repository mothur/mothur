#ifndef FULLMATRIX_H
#define FULLMATRIX_H
/*
 *  fullmatrix.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/6/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "mothur.h"
#include "groupmap.h"
#include "globaldata.hpp"
#include "progress.hpp"


struct Names {
	string		seqName;
	string		groupName;
};

class FullMatrix {
	
public:
	FullMatrix(){};
	FullMatrix(ifstream&);
	~FullMatrix(){};
	
	int getNumSeqs();
	vector<int> getSizes();
	vector<string> getGroups();
	int getNumGroups();
	void printMatrix(ostream&);
	float get(int, int);
	
private:
	vector< vector<float> > matrix;  //a 2D distance matrix of all the sequences and their distances to eachother.
	void readSquareMatrix(ifstream&);  
	void readLTMatrix(ifstream&);
	vector<Names> index; // row in vector, sequence group.  need to know this so when we sort it can be updated.
	vector<int> sizes;
	vector<string> groups;
	void sortGroups(int, int);  //this function sorts the sequences within the matrix.
	
	
	GroupMap* groupmap;  //maps sequences to groups they belong to.
	int numSeqs;
	int numGroups;
	int numUserGroups;
	GlobalData* globaldata;
};

#endif
