#ifndef ORDER_H
#define ORDER_H

#include "datavector.hpp"
#include "sabundvector.hpp"
#include "rabundvector.hpp"

/* This class is a child to datavector.  It represents OTU information at a certain distance. 
	A order vector can be converted into and listvector, rabundvector or sabundvector.
	Each member of the internal container "data" represents the OTU from which it came. 
	So in the example below since there are 6 sequences in OTU 1 there are six 1's in the ordervector.
	and since there are 2 sequences in OTU 3 there are two 3's in the ordervector.
	example: listvector		=	a,b,c,d,e,f		g,h,i		j,k		l		m  
			 rabundvector	=	6				3			2		1		1
			 sabundvector	=	2		1		1		0		0		1
			 ordervector	=	1	1	1	1	1	1	2	2	2	3	3	4	5	*/


class OrderVector : public DataVector {
	
public:
	OrderVector();
//	OrderVector(int);
//	OrderVector(const OrderVector& ov);
	OrderVector(int ns, int nb=0, int mr=0)	: DataVector(), data(ns, -1), maxRank(0), numBins(0), numSeqs(0) {};
	OrderVector(const OrderVector& ov)	: DataVector(ov.label), data(ov.data), maxRank(ov.maxRank), numBins(ov.numBins), numSeqs(ov.numSeqs), needToUpdate(ov.needToUpdate) {if(needToUpdate == 1){	updateStats();}};


	OrderVector(string, vector<int>);
	OrderVector(ifstream&);
	~OrderVector(){};
	
	void set(int, int);
	int get(int);
	void push_back(int);
	void resize(int);
	int size();
	void clear();
	void print(string, ostream&);
	vector<int>::iterator begin();
	vector<int>::iterator end();

	void print(ostream&);

	int getNumBins();
	int getNumSeqs();
	int getMaxRank();
		
	RAbundVector getRAbundVector();
	SAbundVector getSAbundVector();
	OrderVector getOrderVector(map<string,int>*);
	
private:
	vector<int> data;
	int maxRank;
	int numBins;
	int numSeqs;
	bool needToUpdate;
	void updateStats();
};

#endif
