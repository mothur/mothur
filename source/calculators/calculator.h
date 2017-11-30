#ifndef CALCULATOR_H
#define CALCULATOR_H


#include "sabundvector.hpp"
#include "sharedrabundvector.hpp"
#include "sequence.hpp"
#include "mothurout.h"
#include "utils.hpp"

/* The calculator class is the parent class for all the different estimators implemented in mothur except the tree calculators.
It has 2 pure functions EstOutput getValues(SAbundVector*), which works on a single group, and 
EstOutput getValues(SharedRAbundVector* shared1, SharedRAbundVector* shared2), which compares 2 groups. */ 


typedef vector<double> EstOutput;

/***********************************************************************/

class Calculator {

public:
	Calculator(){ m = MothurOut::getInstance(); needsAll = false; }
    Calculator(string n, int c, bool f) : name(n), cols(c), multiple(f) { m = MothurOut::getInstance(); needsAll = false; };
    Calculator(string n, int c, bool f, bool a) : name(n), cols(c), multiple(f), needsAll(a) { m = MothurOut::getInstance(); };
	virtual ~Calculator(){};
	
	virtual EstOutput getValues(SAbundVector*) = 0;	
	virtual EstOutput getValues(vector<SharedRAbundVector*>) = 0;
    //optional calc that returns the otus labels of shared otus
    virtual EstOutput getValues(vector<SharedRAbundVector*> sv , vector<string>&) { data = getValues(sv); return data; }
	virtual void print(ostream& f)	{ f.setf(ios::fixed, ios::floatfield); f.setf(ios::showpoint);
									  f << data[0]; for(int i=1;i<data.size();i++){	f << '\t' << data[i];	}}
    
	virtual string getName()		{	return name;        }
	virtual int getCols()           {	return cols;        }
	virtual bool getMultiple()      {   return multiple;    }
	virtual bool getNeedsAll()      {   return needsAll;    }
    
	virtual string getCitation() = 0;
	void citation() { m->mothurOut(getCitation()); m->mothurOutEndLine(); }
    
protected:
	MothurOut* m;
	EstOutput data;
	string name;
	int cols;
	bool multiple;
	bool needsAll;

};

/**************************************************************************************************/


class ClusterMetric {
    
public:
    ClusterMetric(){ m = MothurOut::getInstance();  }
    ClusterMetric(string n){ m = MothurOut::getInstance();  name = n; }
    virtual ~ClusterMetric(){};
    
    virtual double getValue(long long, long long, long long, long long) = 0; //tp, tn, fp, fn
    
    virtual string getName()		{	return name;        }
    virtual string getCitation() = 0;
    void citation() { m->mothurOut(getCitation()); m->mothurOutEndLine(); }
    
protected:
    MothurOut* m;
    string name;
    
    
};


/**************************************************************************************************/

class DistCalc {
    
public:
    DistCalc(){ dist = 0; m = MothurOut::getInstance(); }
    DistCalc(const DistCalc& d) : dist(d.dist) { m = MothurOut::getInstance(); }
    virtual ~DistCalc() {}
    virtual void calcDist(Sequence, Sequence) = 0;
    double getDist()	{	return dist;	}
    
protected:
    double dist;
    MothurOut* m;
};

/**************************************************************************************************/


#endif

