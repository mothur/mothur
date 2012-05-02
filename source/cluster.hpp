#ifndef CLUSTER_H
#define CLUSTER_H


#include "mothur.h"
#include "sparsematrix.hpp"
#include "mothurout.h"

class RAbundVector;
class ListVector;

typedef vector<MatData> MatVec;

class Cluster {
	
public:
	Cluster(RAbundVector*, ListVector*, SparseMatrix*, float, string);
    virtual void update(double&);				
	virtual string getTag() = 0;
	virtual void setMapWanted(bool m);  
	virtual map<string, int> getSeqtoBin()  {  return seq2Bin;	}

protected:	
	void getRowColCells();
    void removeCell(const MatData& cell, int vrow, int vcol, bool rmMatrix=true);

	virtual bool updateDistance(MatData& colCell, MatData& rowCell) = 0;

	virtual void clusterBins();
	virtual void clusterNames();
	virtual void updateMap();
	
	RAbundVector* rabund;
	ListVector* list;
	SparseMatrix* dMatrix;	
	
	int smallRow;
	int smallCol;
	float smallDist;
	bool mapWanted;
	float cutoff;
	map<string, int> seq2Bin;
	string method;
	
	vector<MatVec> seqVec;		// contains vectors of cells related to a certain sequence
	MatVec rowCells;
	MatVec colCells;
	ull nRowCells;
	ull nColCells;
	MothurOut* m;
};

/***********************************************************************/

class CompleteLinkage : public Cluster {
public:
	CompleteLinkage(RAbundVector*, ListVector*, SparseMatrix*, float, string);
	bool updateDistance(MatData& colCell, MatData& rowCell);
	string getTag();
	
private:
		
};

/***********************************************************************/

class SingleLinkage : public Cluster {
public:
	SingleLinkage(RAbundVector*, ListVector*, SparseMatrix*, float, string);
    void update(double&);
	bool updateDistance(MatData& colCell, MatData& rowCell);
	string getTag();
	
private:
		
};

/***********************************************************************/

class AverageLinkage : public Cluster {
public:
	AverageLinkage(RAbundVector*, ListVector*, SparseMatrix*, float, string);
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

class WeightedLinkage : public Cluster {
public:
	WeightedLinkage(RAbundVector*, ListVector*, SparseMatrix*, float, string);
	bool updateDistance(MatData& colCell, MatData& rowCell);
	string getTag();
	
private:
	int saveRow;
	int saveCol;	
};

/***********************************************************************/

#endif
