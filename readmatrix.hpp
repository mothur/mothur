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

#include "rabundvector.hpp"
#include "listvector.hpp"
#include "sparsematrix.hpp"
#include "nameassignment.hpp"
#include "inputdata.h"
#include "globaldata.hpp"
#include "sabundvector.hpp"
#include "groupmap.h"

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
	OrderVector* order;
	InputData* input;
	float cutoff;
};



class ReadPhylipMatrix : public ReadMatrix {
	
public:
	ReadPhylipMatrix(string);
	~ReadPhylipMatrix();
	void read(NameAssignment*);
private:
	ifstream fileHandle;
	string distFile;
};



class ReadColumnMatrix : public ReadMatrix {
	
public:
	ReadColumnMatrix(string);
	~ReadColumnMatrix();
	void read(NameAssignment*);
private:
	ifstream fileHandle;
	string distFile;
};


class ReadPhilFile : public ReadMatrix {
	
public:
	ReadPhilFile(string);
	~ReadPhilFile();
	void read(GlobalData* globaldata);
private:
	ifstream fileHandle;
	string philFile;
	InputData* input;
	InputData* inputSabund;
	ListVector* list;
	OrderVector* order;
	SAbundVector* sabund;
	GlobalData* globaldata;
	// InputData* getInput()			{	return input;	}
};



#endif
