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

class FullMatrix {
	
	public:
		FullMatrix(){};
		FullMatrix(ifstream&);
		~FullMatrix(){};
	
		int getNumSeqs();
		void printMatrix(ostream&);
	
	private:
		void sortGroups();  //this function sorts the sequences within the matrix.
		void quicksort(int, int, int);//row of matrix, low, high and row number
		void readSquareMatrix(ifstream&);  
		void readLTMatrix(ifstream&);
		vector< vector<float> > matrix;  //a 2D distance matrix of all the sequences and their distances to eachother.
		map<int, string> index; // row in vector, sequence group.  need to know this so when we sort it can be updated.
		GroupMap* groupmap;  //maps sequences to groups they belong to.
		GlobalData* globaldata;
		int numSeqs;
		bool square;

};

#endif