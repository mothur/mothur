#ifndef SABUND_H
#define SABUND_H

using namespace std;

#include "datavector.hpp"
#include "rabundvector.hpp"
#include "ordervector.hpp"

class SAbundVector : public DataVector {
	
public:
	SAbundVector();
	SAbundVector(int);
//	SAbundVector(const SAbundVector&);
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
	void resize(int);
	int size();

	void print(ostream&);
	void print(string, ostream&);
		
	RAbundVector getRAbundVector();	
	SAbundVector getSAbundVector();
	OrderVector getOrderVector(map<string,int>*);
	
private:
	vector<int> data;
//	bool needToUpdate;
//	void updateStats();
	int maxRank;
	int numBins;
	int numSeqs;

};

#endif
