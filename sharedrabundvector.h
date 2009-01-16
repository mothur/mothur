#ifndef SHAREDVECTOR_H
#define SHAREDVECTOR_H

/*
 *  sharedvector.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 12/5/08.
 *  Copyright 2008 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include <Carbon/Carbon.h>
#include "datavector.hpp"
#include "sharedordervector.h"




class SharedRAbundVector : public DataVector {
	
public:
	SharedRAbundVector();
	SharedRAbundVector(int);
	//SharedRAbundVector(string, vector<int>);
	SharedRAbundVector(const SharedRAbundVector& bv) : DataVector(bv), data(bv.data), maxRank(bv.maxRank), numBins(bv.numBins), numSeqs(bv.numSeqs){};
//	SharedRAbundVector(ifstream&);
	~SharedRAbundVector();

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
	vector<individual>::reverse_iterator rbegin();
	vector<individual>::reverse_iterator rend();
	
	void print(ostream&);
		
	SharedRAbundVector getSharedRAbundVector();
	SAbundVector getSAbundVector();
	OrderVector getOrderVector(map<string,int>*);
	SharedOrderVector getSharedOrderVector();
	
private:
	vector<individual>  data; 
	int maxRank;
	int numBins;
	int numSeqs;
	string group;	
};


#endif

