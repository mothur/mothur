#ifndef CLUSTER_H
#define CLUSTER_H


#include "mothur.h"
#include "sparsematrix.hpp"

class RAbundVector;
class ListVector;

typedef vector<MatData> MatVec;

class Cluster {
	
public:
	Cluster(RAbundVector*, ListVector*, SparseMatrix*);
    virtual void update();
	virtual string getTag() = 0;

protected:	
	void getRowColCells();
    void removeCell(const MatData& cell, int vrow, int vcol, bool rmMatrix=true);

	virtual bool updateDistance(MatData& colCell, MatData& rowCell) = 0;

	virtual void clusterBins();
	virtual void clusterNames();
	
	RAbundVector* rabund;
	ListVector* list;
	SparseMatrix* dMatrix;	
	
	int smallRow;
	int smallCol;
	float smallDist;
	
	vector<MatVec> seqVec;		// contains vectors of cells related to a certain sequence
	MatVec rowCells;
	MatVec colCells;
	ull nRowCells;
	ull nColCells;
};

/***********************************************************************/

class CompleteLinkage : public Cluster {
public:
	CompleteLinkage(RAbundVector*, ListVector*, SparseMatrix*);
	bool updateDistance(MatData& colCell, MatData& rowCell);
	string getTag();
	
private:
		
};

/***********************************************************************/

class SingleLinkage : public Cluster {
public:
	SingleLinkage(RAbundVector*, ListVector*, SparseMatrix*);
    void update();
	bool updateDistance(MatData& colCell, MatData& rowCell);
	string getTag();
	
private:
		
};

/***********************************************************************/

class AverageLinkage : public Cluster {
public:
	AverageLinkage(RAbundVector*, ListVector*, SparseMatrix*);
	bool updateDistance(MatData& colCell, MatData& rowCell);
	string getTag();
	
private:
	int saveRow;
	int saveCol;
	int rowBin;
	int colBin;
	int totalBin;

};

/***********************************************************************/

#endif
