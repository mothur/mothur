#ifndef SABUND_H
#define SABUND_H

#include "datavector.hpp"
#include "rabundvector.hpp"
#include "ordervector.hpp"
#include "calculator.h"


/*  Data Structure for a sabund file.
    This class is a child to datavector.  It represents OTU information at a certain distance. 
	A sabundvector can be converted into and ordervector, listvector or rabundvector.
	Each member of the internal container "data" represents the number of OTU's with that many members, but staring at 1.
	So data[1] = 2, because there are two OTUs with 1 member.
	example: listvector		=	a,b,c,d,e,f		g,h,i		j,k		l		m  
			 rabundvector	=	6				3			2		1		1
			 sabundvector	=	2		1		1		0		0		1
			 ordervector	=	1	1	1	1	1	1	2	2	2	3	3	4	5 */


class SAbundVector : public DataVector {
	
public:
	SAbundVector();
	SAbundVector(int);
//	SAbundVector(const SAbundVector&);
	SAbundVector(vector<int>, int, int, int);
	SAbundVector(string, vector<int>);
	SAbundVector(const SAbundVector& rv) : DataVector(rv.label), data(rv.data), maxRank(rv.maxRank), numBins(rv.numBins), numSeqs(rv.numSeqs){};
	SAbundVector(ifstream&);
	~SAbundVector(){};

	int getNumBins();	
	int getNumSeqs();	
	int getMaxRank();	
	
	void set(int, int);
	int get(int);
	void push_back(int);
	void quicksort();
	int sum();
	void resize(int);
	int size();
	void clear();

	void print(ostream&);
	void print(string, ostream&);
		
	RAbundVector getRAbundVector();		
	SAbundVector getSAbundVector();
	OrderVector getOrderVector(map<string,int>* hold = nullptr);
	
private:
	vector<int> data;
//	bool needToUpdate;
//	void updateStats();
	
	int maxRank;
	int numBins;
	int numSeqs;

};

#endif
