#ifndef CLUSTER_H
#define CLUSTER_H

#include "sparsedistancematrix.h"
#include "optimatrix.h"
#include "mothurout.h"
#include "rabundvector.hpp"
#include "listvector.hpp"

class ListVector;

class Cluster {
	
public:
	Cluster(RAbundVector*, ListVector*, SparseDistanceMatrix*, float, string, float);
    Cluster(){}
    virtual ~Cluster() {}
    virtual bool update(double&);
	virtual string getTag() = 0;
	virtual void setMapWanted(bool m);  
	virtual map<string, int> getSeqtoBin()  {  return seq2Bin;	}
    
protected:	    
	virtual bool updateDistance(PDistCell& colCell, PDistCell& rowCell) = 0;
    
	virtual void clusterBins();
	virtual void clusterNames();
	virtual void updateMap();
	
	RAbundVector* rabund;
	ListVector* list;
	SparseDistanceMatrix* dMatrix;	
	
	ull smallRow;
	ull smallCol;
	float smallDist, adjust;
	bool mapWanted;
	float cutoff;
	map<string, int> seq2Bin;
	string method;
	
	ull nRowCells;
	ull nColCells;
	MothurOut* m;
};

/***********************************************************************/

class CompleteLinkage : public Cluster {
public:
	CompleteLinkage(RAbundVector*, ListVector*, SparseDistanceMatrix*, float, string, float);
	bool updateDistance(PDistCell& colCell, PDistCell& rowCell);
	string getTag();
	
private:
    
};

/***********************************************************************/

class SingleLinkage : public Cluster {
public:
	SingleLinkage(RAbundVector*, ListVector*, SparseDistanceMatrix*, float, string, float);
    //void update(double&);
	bool updateDistance(PDistCell& colCell, PDistCell& rowCell);
	string getTag();
	
private:
    
};

/***********************************************************************/

class AverageLinkage : public Cluster {
public:
	AverageLinkage(RAbundVector*, ListVector*, SparseDistanceMatrix*, float, string, float);
	bool updateDistance(PDistCell& colCell, PDistCell& rowCell);
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
	WeightedLinkage(RAbundVector*, ListVector*, SparseDistanceMatrix*, float, string, float);
	bool updateDistance(PDistCell& colCell, PDistCell& rowCell);
	string getTag();
	
private:
	int saveRow;
	int saveCol;	
};

/***********************************************************************/



#endif
