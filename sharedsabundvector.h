#ifndef SHAREDSABUND_H
#define SHAREDSABUND_H


/*
 *  sharedSharedSAbundVector.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 12/10/08.
 *  Copyright 2008 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "datavector.hpp"
#include "rabundvector.hpp"
#include "ordervector.hpp"
#include "sharedordervector.h"
#include "sharedrabundvector.h"

/* This class is a child to datavector.  It represents OTU information at a certain distance. 
	It is similiar to an sabundvector except each member of data knows which group it belongs to.
	Each member of the internal container "data" is a struct of type individual. 
	An individual which knows the OTU from which it came, 
	the group it is in and its abundance.  */


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
	
	void set(int, int, string);	 //OTU, abundance, group
	individual get(int);
	int getAbundance(int);
	void push_back(int, int, string);	//abundance, OTU, group
	void pop_back();
	void resize(int);
	int size();
	void clear();

	void print(ostream&);
		
	RAbundVector getRAbundVector();	
	SAbundVector getSAbundVector();
	OrderVector getOrderVector(map<string,int>*);
	SharedSAbundVector getSharedSAbundVector();
	SharedRAbundVector getSharedRAbundVector();
	SharedOrderVector getSharedOrderVector();
	
private:
	vector<individual> data;
	
	int maxRank;
	int numBins;
	int numSeqs;
	string group;
};

#endif

