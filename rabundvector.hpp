#ifndef RABUND_H
#define RABUND_H

#include "datavector.hpp"

class RAbundVector : public DataVector {
	
public:
	RAbundVector();
	RAbundVector(int);
//	RAbundVector(const RAbundVector&);
	RAbundVector(string, vector<int>);
	RAbundVector(const RAbundVector& bv) : DataVector(bv), data(bv.data), maxRank(bv.maxRank), numBins(bv.numBins), numSeqs(bv.numSeqs){};
	RAbundVector(ifstream&);
	~RAbundVector();

	int getNumBins();		
	int getNumSeqs();							
	int getMaxRank();							

	void set(int, int);	
	int get(int);
	void push_back(int);
	void pop_back();
	void resize(int);
	int size();
	vector<int>::reverse_iterator rbegin();
	vector<int>::reverse_iterator rend();
	
	void print(ostream&);
	void print(string, ostream&);
	
	RAbundVector getRAbundVector();
	SAbundVector getSAbundVector();
	OrderVector getOrderVector(map<string,int>*);
	
private:
	vector<int> data;
	int maxRank;
	int numBins;
	int numSeqs;	
};


#endif
