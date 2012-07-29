#ifndef READMATRIX_HPP
#define READMATRIX_HPP

/*
 *  readmatrix.hpp
 *  
 *
 *  Created by Pat Schloss on 8/13/08.
 *  Copyright 2008 Patrick D. Schloss. All rights reserved.
 *
 */

#include "mothur.h"
#include "listvector.hpp"
#include "nameassignment.hpp"
#include "sparsedistancematrix.h"

class SparseMatrix;

class ReadMatrix {

public:
	ReadMatrix(){ DMatrix = new SparseDistanceMatrix(); m = MothurOut::getInstance();  }
	virtual ~ReadMatrix() {}
	virtual int read(NameAssignment*){ return 1; }
	
	void setCutoff(float c)			{	cutoff = c;		}
    SparseDistanceMatrix* getDMatrix()		{	return DMatrix;		}
	ListVector* getListVector()		{	return list;	}

	int successOpen;
	
protected:
    SparseDistanceMatrix* DMatrix;
	ListVector* list;
	float cutoff;
	MothurOut* m;
	bool sim;
};



#endif
