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

class FullMatrix {
	
public:
	FullMatrix(){};
	FullMatrix(ifstream&);
	~FullMatrix();
	
	int getNumSeqs();
	
private:
	void sortGroups();  //this function sorts the sequences within the matrix.  
	vector< vector<float> > matrix;  //a 2D distance matrix of all the sequences and their distances to eachother.
	GroupMap* groupmap;  //maps sequences to groups they belong to.
	int numSeqs;

};

#endif