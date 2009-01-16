#ifndef ORDER_H
#define ORDER_H

#include "datavector.hpp"
#include "sabundvector.hpp"
#include "rabundvector.hpp"


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
