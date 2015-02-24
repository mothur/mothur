#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include "mothur.h"
#include "mothurout.h"


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

typedef list<PCell>::iterator MatData;

class SparseMatrix {
	
public:
	SparseMatrix();
	~SparseMatrix(){  while(!mins.empty() && mins.back() == NULL){  mins.pop_back();	}  }
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
};

/***********************************************************************/

#endif
