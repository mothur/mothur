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
#include "progress.hpp"


struct Names {
	string		seqName;
	string		groupName;
};

class FullMatrix {
	
public:
	//FullMatrix(){ m = MothurOut::getInstance(); }
	FullMatrix(ifstream&, GroupMap*, bool);
	~FullMatrix(){};
	
	int getNumSeqs();
	vector<int> getSizes();
	vector<string> getGroups();
	void setGroups(vector<string> names) { groups = names;  }
	void setSizes(vector<int> s)		 { sizes = s;		}
	int getNumGroups();
	void printMatrix(ostream&);
	float get(int, int);
	Names getRowInfo(int row)  {  return index[row];  }
	
private:
	vector< vector<float> > matrix;  //a 2D distance matrix of all the sequences and their distances to eachother.
	int readSquareMatrix(ifstream&);  
	int readLTMatrix(ifstream&);
	vector<Names> index; // row in vector, sequence group.  need to know this so when we sort it can be updated.
	vector<int> sizes;
	vector<string> groups;
	
	void sortGroups(int, int);  //this function sorts the sequences within the matrix.
	void swapRows(int, int);
	
	GroupMap* groupmap;  //maps sequences to groups they belong to.
	int numSeqs;
	int numGroups;
	int numUserGroups;
	bool sim;
	MothurOut* m;
};

#endif
