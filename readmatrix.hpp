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
#include "sparsematrix.hpp"
#include "nameassignment.hpp"
#include "globaldata.hpp"

class SparseMatrix;

class ReadMatrix {

public:
	ReadMatrix(){	D = new SparseMatrix();	}
	virtual void read(NameAssignment*){};
	virtual void read(GlobalData* globaldata){};
	void setCutoff(float c)			{	cutoff = c;		}
	SparseMatrix* getMatrix()		{	return D;		}
	ListVector* getListVector()		{	return list;	}
//	OrderVector* getOrderVector()   {	return order;	}

	int successOpen;
	
protected:
	SparseMatrix* D;
	ListVector* list;
	GlobalData* globaldata;
	float cutoff;
};



#endif
