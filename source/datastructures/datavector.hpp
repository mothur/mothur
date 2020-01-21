#ifndef datavector_h
#define datavector_h

#include "mothurout.h"
#include "utils.hpp"

/* This class is parent to listvector, ordervector, rabundvector, sabundvector, sharedordervector, sharedrabundvector, sharedsabundvector. 
	The child classes all contain OTU information in different forms. */
	
class SharedSAbundVector;
class SharedRAbundVectors;
class SharedLCRVector;
class SharedLCRVectors;
class RAbundVector;
class RAbundFloatVector;
class SAbundVector;
class OrderVector;
class SharedListVector;
class SharedOrderVector;
class GroupMap;

class DataVector {
	
public:
	DataVector(){ m = MothurOut::getInstance(); }
	DataVector(string l) : label(l) { m = MothurOut::getInstance();};
    DataVector(const DataVector& dv) : label(dv.label){ m = MothurOut::getInstance();}
	DataVector(ifstream&) {m = MothurOut::getInstance();}
	DataVector(ifstream&, GroupMap*){m = MothurOut::getInstance();}
	virtual ~DataVector(){};

	virtual int size() = 0;
    virtual void clear() = 0;
	
    virtual RAbundVector getRAbundVector() = 0;
	virtual SAbundVector getSAbundVector() = 0;
    virtual OrderVector getOrderVector(map<string,int>* hold = NULL) = 0;
    virtual void resize(int) = 0;
    
    virtual void print(ostream&, map<string, int>&) {}
    virtual void print(ostream&, bool) { m->mothurOut("[ERROR]: no print function\n"); }
    virtual void printHeaders(ostream&) {};
    virtual void print(ostream&, bool&) {}
    virtual void print(ostream&) {}
    
    void setLabel(string l)		{	label = l;			}
    string getLabel()			{	return label;		}
    
	
protected:
	string label;
	MothurOut* m;
    Utils util;

};

/***********************************************************************/

#endif
