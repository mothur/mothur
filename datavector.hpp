#ifndef datavector_h
#define datavector_h


#include "mothur.h"
#include "mothurout.h"

/* This class is parent to listvector, ordervector, rabundvector, sabundvector, sharedordervector, sharedrabundvector, sharedsabundvector. 
	The child classes all contain OTU information in different forms. */
	

class RAbundVector;
class SAbundVector;
class OrderVector;
class SharedListVector;
class SharedOrderVector;
class SharedSAbundVector;
class SharedRAbundVector;
class SharedRAbundFloatVector;
class GroupMap;

class DataVector {
	
public:
	DataVector(){ m = MothurOut::getInstance(); }// : maxRank(0), numBins(0), numSeqs(0){};
	DataVector(string l) : label(l) { m = MothurOut::getInstance();};
	DataVector(const DataVector& dv) : label(dv.label){ m = MothurOut::getInstance();};//, maxRank(dv.maxRank), numBins(dv.numBins), numSeqs(dv.numSeqs) {};
	DataVector(ifstream&) {m = MothurOut::getInstance();}
	DataVector(ifstream&, GroupMap*){m = MothurOut::getInstance();}
	virtual ~DataVector(){};
	
//	virtual int getNumBins()	{	return numBins;		}
//	virtual int getNumSeqs()	{	return numSeqs;		}
//	virtual int getMaxRank()	{	return maxRank;		}
	
	virtual void resize(int) = 0;
	virtual int size()	= 0;
	virtual void print(ostream&);
    virtual void print(ostream&, map<string, int>&) {}
	virtual void printHeaders(ostream&) {};
	virtual void clear() = 0;
	
	void setLabel(string l)		{	label = l;			}
	string getLabel()			{	return label;		}

	virtual RAbundVector getRAbundVector() = 0;
	virtual SAbundVector getSAbundVector() = 0;
	virtual OrderVector getOrderVector(map<string,int>* hold = NULL) = 0;
	
protected:
	string label;
	MothurOut* m;
//	int maxRank;
//	int numBins;
//	int numSeqs;	
};

/***********************************************************************/

#endif
