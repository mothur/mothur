#ifndef LIST_H
#define LIST_H

#include "datavector.hpp"
#include "groupmap.h"
#include "globaldata.hpp"
#include "sharedordervector.h"
#include <iostream>
#include <map>


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
	SharedOrderVector* getSharedOrderVector();
	
private:
	vector<string> data;
	GlobalData* globaldata;
	GroupMap* groupmap;
	int maxRank;
	int numBins;
	int numSeqs;

};

#endif
