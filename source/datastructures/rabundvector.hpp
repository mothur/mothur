#ifndef RABUND_H
#define RABUND_H

#include "datavector.hpp"

/*  Data Structure for a rabund file.
    This class is a child to datavector.  It represents OTU information at a certain distance. 
	A rabundvector can be converted into and ordervector, listvector or sabundvector.
	Each member of the internal container "data" represents an individual OTU.
	So data[0] = 6, because there are six member in that OTU.
	example: listvector		=	a,b,c,d,e,f		g,h,i		j,k		l		m  
			 rabundvector	=	6				3			2		1		1
			 sabundvector	=	2		1		1		0		0		1
			 ordervector	=	1	1	1	1	1	1	2	2	2	3	3	4	5 */

class RAbundFloatVector;
//class OrderVector;

class RAbundVector : public DataVector {
	
public:
	RAbundVector();
	RAbundVector(int);
	RAbundVector(vector<int>, int, int, int);
	RAbundVector(string, vector<int>);
	RAbundVector(const RAbundVector& bv) : DataVector(bv), data(bv.data), maxRank(bv.maxRank), numBins(bv.numBins), numSeqs(bv.numSeqs){};
	RAbundVector(ifstream&);
    RAbundVector(ifstream& f, string l); //label given
	~RAbundVector();

	int getNumBins();		
	int getNumSeqs();							
	int getMaxRank();							

    int remove(int);
	void set(int, int);	
	int get(int);
	void push_back(int);
	void pop_back();
	void resize(int);
	int size();
	void quicksort();
	int sum();
	int sum(int);
	int numNZ();
    vector<int> getSortedD();
	void clear();
	vector<int>::reverse_iterator rbegin();
	vector<int>::reverse_iterator rend();
	
	void print(ostream&); //sorted, no group
    int print(ostream&, string); //notsorted, group
	void print(string, ostream&); //label, sorted
	void nonSortedPrint(ostream&); //nonsorted , no group
	
	RAbundVector getRAbundVector();
    RAbundFloatVector getRAbundFloatVector();
	SAbundVector getSAbundVector();
	OrderVector getOrderVector(map<string,int>*);
	
private:
	vector<int> data;
	int maxRank;
	int numBins;
	int numSeqs;	
};


#endif
