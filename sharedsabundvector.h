#ifndef SHAREDSABUND_H
#define SHAREDSABUND_H


/*
 *  sharedSharedSAbundVector.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 12/10/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */


#include <Carbon/Carbon.h>
#include "datavector.hpp"
#include "rabundvector.hpp"
#include "ordervector.hpp"
#include "sharedordervector.h"
#include "sharedrabundvector.h"

using namespace std;

class SharedSAbundVector : public DataVector {
	
public:
	SharedSAbundVector();
	SharedSAbundVector(int);
	SharedSAbundVector(const SharedSAbundVector& rv) : DataVector(rv.label), data(rv.data), maxRank(rv.maxRank), numBins(rv.numBins), numSeqs(rv.numSeqs){};
	~SharedSAbundVector(){};

	int getNumBins();	
	int getNumSeqs();	
	int getMaxRank();	
	string getGroup();
	void setGroup(string);	
	
	void set(int, int, string);	
	individual get(int);
	int getAbundance(int);
	void push_back(int, int, string);
	void pop_back();
	void resize(int);
	int size();

	void print(ostream&);
		
	RAbundVector getRAbundVector();	
	SAbundVector getSAbundVector();
	SharedSAbundVector getSharedSAbundVector();
	SharedRAbundVector getSharedVector();
	OrderVector getOrderVector(map<string,int>*);
	
private:
	vector<individual> data;
	
	int maxRank;
	int numBins;
	int numSeqs;
	string group;
};

#endif

