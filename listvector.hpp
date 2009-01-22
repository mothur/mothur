#ifndef LIST_H
#define LIST_H

#include <Carbon/Carbon.h>
#include "datavector.hpp"
#include <iostream>
#include <map>

/* This class is a child to datavector.  It represents OTU information at a certain distance. 
	A list vector can be converted into and ordervector, rabundvector or sabundvector.
	Each member of the internal container "data" represents an individual OTU.
	So data[0] = "a,b,c,d,e,f".
	example: listvector		=	a,b,c,d,e,f		g,h,i		j,k		l		m  
			 rabundvector	=	6				3			2		1		1
			 sabundvector	=	2		1		1		0		0		1
			 ordervector	=	1	1	1	1	1	1	2	2	2	3	3	4	5 */

class ListVector : public DataVector {
	
public:
	ListVector();
	ListVector(int);
//	ListVector(const ListVector&);
	ListVector(string, vector<string>);
	ListVector(const ListVector& lv) : DataVector(lv.label), data(lv.data), maxRank(lv.maxRank), numBins(lv.numBins), numSeqs(lv.numSeqs){};
	ListVector(ifstream&);
	~ListVector(){};
	
	int getNumBins()							{	return numBins;		}
	int getNumSeqs()							{	return numSeqs;		}
	int getMaxRank()							{	return maxRank;		}

	void set(int, string);	
	string get(int);
	void push_back(string);
	void resize(int);
	void clear();
	int size();
	void print(ostream&);
	
	RAbundVector getRAbundVector();
	SAbundVector getSAbundVector();
	OrderVector getOrderVector(map<string,int>*);
	
private:
	vector<string> data;  //data[i] is a list of names of sequences in the ith OTU.
	int maxRank;
	int numBins;
	int numSeqs;

};

#endif
