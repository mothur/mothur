#ifndef CLUSTER_H
#define CLUSTER_H

using namespace std;

#include <iostream>
#include <list>

#include "utilities.hpp"
#include "sparsematrix.hpp"
#include "rabundvector.hpp"

class RAbundVector;
class ListVector;
class SparseMatrix;

typedef list<PCell>::iterator MatData;

class Cluster {
	
public:
	Cluster(RAbundVector*, ListVector*, SparseMatrix*);
	virtual void update() = 0;
	
protected:	
	void getRowColCells();
	virtual void clusterBins();
	virtual void clusterNames();
	
	RAbundVector* rabund;
	ListVector* list;
	SparseMatrix* dMatrix;	
	
	int smallRow;
	int smallCol;
	float smallDist;
	vector<MatData> rowCells;
	vector<MatData> colCells;
	ull nRowCells;
	ull nColCells;
};

/***********************************************************************/

class CompleteLinkage : public Cluster {
public:
	CompleteLinkage(RAbundVector*, ListVector*, SparseMatrix*);
	void update();
	
private:
		
};

/***********************************************************************/

class SingleLinkage : public Cluster {
public:
	SingleLinkage(RAbundVector*, ListVector*, SparseMatrix*);
	void update();
	
private:
		
};

/***********************************************************************/

class AverageLinkage : public Cluster {
public:
	AverageLinkage(RAbundVector*, ListVector*, SparseMatrix*);
	void update();
	
private:
		
};

/***********************************************************************/

#endif
