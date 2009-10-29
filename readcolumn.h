#ifndef READCOLUMN_H
#define READCOLUMN_H
/*
 *  readcolumn.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/21/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "readmatrix.hpp"

/******************************************************/

class ReadColumnMatrix : public ReadMatrix {
	
public:
	ReadColumnMatrix(string);
	~ReadColumnMatrix();
	void read(NameAssignment*);
private:
	ifstream fileHandle;
	string distFile;
	
};

/******************************************************/

#endif
