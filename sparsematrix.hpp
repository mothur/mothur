#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include "mothur.h"


class ListVector;

/***********************************************************************/


 struct PCell{
	ull row;
	ull column;
	float dist;
	PCell** vectorMap;
	PCell() : row(0), column(0), dist(0), vectorMap(NULL) {};
	PCell(ull r, ull c, float d) : row(r), column(c), dist(d), vectorMap(NULL) {};
};

/***********************************************************************/

class SparseMatrix {
	
public:
	SparseMatrix();
	int getNNodes();
	void print();					//Print the contents of the matrix
	void print(ListVector*);					//Print the contents of the matrix
	PCell* getSmallestCell();		//Return the cell with the smallest distance
	float getSmallDist();
	
	void rmCell(list<PCell>::iterator);
	void addCell(PCell);
	void clear();
	list<PCell>::iterator begin();
	list<PCell>::iterator end();

private:
	PCell* smallCell;				//The cell with the smallest distance
	int numNodes;

	list<PCell> matrix;
	
	vector<PCell*> mins;
	float smallDist;
	int minsIndex;
};

/***********************************************************************/

#endif
