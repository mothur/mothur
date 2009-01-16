#ifndef SHAREDORDER_H
#define SHAREDORDER_H
/*
 *  sharedorder.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 12/9/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <Carbon/Carbon.h>
#include "datavector.hpp"
#include "sabundvector.hpp"
#include "rabundvector.hpp"

struct individual {
		string group;
		int bin;
		int abundance;
};


class SharedOrderVector : public DataVector {
	
public:
	SharedOrderVector();
//	SharedOrderVector(int ns, int nb=0, int mr=0)	: DataVector(), data(ns, -1), maxRank(0), numBins(0), numSeqs(0) {};
	SharedOrderVector(const SharedOrderVector& ov)	: DataVector(ov.label), data(ov.data), maxRank(ov.maxRank), numBins(ov.numBins), numSeqs(ov.numSeqs), needToUpdate(ov.needToUpdate) {if(needToUpdate == 1){	updateStats();}};

	SharedOrderVector(string, vector<individual>);
//	SharedOrderVector(ifstream&);
	~SharedOrderVector(){};
	
	void set(int, int, int, string);
	individual get(int);
	void push_back(int, int, string);
	void resize(int);
	int size();
	void print(ostream&);
	vector<individual>::iterator begin();
	vector<individual>::iterator end();


	int getNumBins();
	int getNumSeqs();
	int getMaxRank();
		
	RAbundVector getRAbundVector();
	SAbundVector getSAbundVector();
	OrderVector getOrderVector(map<string,int>*);
	SharedOrderVector getSharedOrderVector();
	
private:
	vector<individual>  data; 
	map< int, vector<individual> >::iterator it;
	int maxRank;
	int numBins;
	int numSeqs;
	bool needToUpdate;
	void updateStats();
};

#endif

