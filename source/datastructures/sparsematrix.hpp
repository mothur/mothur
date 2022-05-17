#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include "mothur.h"
#include "mothurout.h"
#include "utils.hpp"


class ListVector;

/***********************************************************************/

typedef list<PCell>::iterator MatData;

class SparseMatrix {
	
public:
	SparseMatrix();
	~SparseMatrix(){  while(!mins.empty() && mins.back() == nullptr){  mins.pop_back();	}  }
	int getNNodes();
	void print();					//Print the contents of the matrix
	void print(ListVector*);		//Print the contents of the matrix
	PCell* getSmallestCell();		//Return the cell with the smallest distance
	float getSmallDist();
	
	MatData rmCell(MatData);
	void addCell(PCell);
	void clear();
	MatData begin();
	MatData end();

private:
	PCell* smallCell;				//The cell with the smallest distance
	int numNodes;

	list<PCell> matrix;
	
	vector<PCell*> mins;
	float smallDist;
	int minsIndex;
	MothurOut* m;
    Utils util;
};

/***********************************************************************/

#endif
