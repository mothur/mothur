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

struct Swap {
	int		a;
	int		b;
};


class FullMatrix {
	
	public:
		FullMatrix(){};
		FullMatrix(ifstream&);
		~FullMatrix(){};
	
		int getNumSeqs();
		void printMatrix(ostream&);
		void setBounds();  //requires globaldata->Groups to be filled
		vector<float> getMins(int); //returns vector of mins for "box" requested ie. groups A, B, 0 = AA, 1 = AB, 2 = BA, 3 = BB;
		void getDist(vector<float>&);  //fills a vector with the valid distances for the integral form.
		void shuffle(string, string);  //shuffles the sequences in the groups passed in.
		void restore();  //unshuffles the matrix.
	
	private:
		void sortGroups(int, int);  //this function sorts the sequences within the matrix.
		void getBounds(int&, string);
		void readSquareMatrix(ifstream&);  
		void readLTMatrix(ifstream&);
		void printMinsForRows(ostream&);
		
		map<int, Names> index; // row in vector, sequence group.  need to know this so when we sort it can be updated.
		map<int, Swap> restoreIndex; //a map of the swaps made so you can undo them in restore.
		map<int, Names>::iterator it;
		map<int, Swap>::reverse_iterator it2;
			
		vector< vector<float> > matrix;  //a 2D distance matrix of all the sequences and their distances to eachother.
		vector<float> minsForRows;  //vector< minimum distance for that subrow> - one for each comparison.
		vector<int> bounds;  //bounds[1] = starting row in matrix from group B, bounds[2] = starting row in matrix from group C, bounds[3] = no need to find upper bound of C because its numSeqs.
								
		GroupMap* groupmap;  //maps sequences to groups they belong to.
		GlobalData* globaldata;
		int numSeqs, numGroups, numUserGroups;
		bool square;

};

#endif