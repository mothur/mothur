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

using namespace std;

struct Names {
	string		groupname;
	string		seqName;
};


class FullMatrix {
	
	public:
		FullMatrix(){};
		FullMatrix(ifstream&);
		~FullMatrix(){};
	
		int getNumSeqs();
		void printMatrix(ostream&);
		void getMinsForRowsVectors();  //requires globaldata->Groups to be filled
	
	private:
		void sortGroups(int, int);  //this function sorts the sequences within the matrix.
		void getBounds(int&, string);
		void readSquareMatrix(ifstream&);  
		void readLTMatrix(ifstream&);
		vector< vector<float> > matrix;  //a 2D distance matrix of all the sequences and their distances to eachother.
		vector< vector<float> > minsForRows;  //vector< minimum distance for that subrow> -one for each comparison.
		map<int, Names> index; // row in vector, sequence group.  need to know this so when we sort it can be updated.
		map<int, Names>::iterator it;
		GroupMap* groupmap;  //maps sequences to groups they belong to.
		GlobalData* globaldata;
		int numSeqs, numGroups, numUserGroups;
		bool square;

};

#endif