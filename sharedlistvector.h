#ifndef SHAREDLIST_H
#define SHAREDLIST_H

/*
 *  sharedlistvector.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/22/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "datavector.hpp"
#include "groupmap.h"
#include "counttable.h"
#include "sharedrabundvector.h"
#include "sharedsabundvector.h"

/* This class is a child to datavector.  It represents OTU information at a certain distance. 
	A sharedlistvector can be converted into a sharedordervector, sharedrabundvector or sharedsabundvectorand 
	as well as an ordervector, rabundvector or sabundvector.
	Each member of the internal container "data" represents an individual OTU.
	Each individual in the OTU belongs to a group.
	So data[0] = "a,b,c,d,e,f".
	example: listvector		=	a,b,c,d,e,f		g,h,i		j,k		l		m  
			 rabundvector	=	6				3			2		1		1
			 sabundvector	=	2		1		1		0		0		1
			 ordervector	=	1	1	1	1	1	1	2	2	2	3	3	4	5 */

class SharedListVector : public DataVector {
	
public:
	SharedListVector();
	SharedListVector(int);
	SharedListVector(ifstream&);
	SharedListVector(const SharedListVector& lv) : DataVector(lv.label), data(lv.data), maxRank(lv.maxRank), numBins(lv.numBins), numSeqs(lv.numSeqs), binLabels(lv.binLabels) { groupmap = NULL; countTable = NULL; };
	~SharedListVector(){ if (groupmap != NULL) { delete groupmap; } if (countTable != NULL) { delete countTable; } };
	
	int getNumBins()							{	return numBins;		}
	int getNumSeqs()							{	return numSeqs;		}
	int getMaxRank()							{	return maxRank;		}

	void set(int, string);	
	string get(int);
    vector<string> getLabels();
    void setLabels(vector<string>);
	void push_back(string);
	void resize(int);
	void clear();
	int size();
	void print(ostream&);
	
	RAbundVector getRAbundVector();
	SAbundVector getSAbundVector();
	OrderVector getOrderVector(map<string,int>*);
	SharedOrderVector* getSharedOrderVector();
	SharedRAbundVector getSharedRAbundVector(string);  //get sharedrabundvector for a certain group
	SharedSAbundVector getSharedSAbundVector(string);			//get sharedsabundvector for a certain group
	vector<SharedRAbundVector*> getSharedRAbundVector(); //returns sharedRabundVectors for all the users groups
	
private:
	vector<string> data;  //data[i] is a list of names of sequences in the ith OTU.
	GroupMap* groupmap;
    CountTable* countTable;
	int maxRank;
	int numBins;
	int numSeqs;
    vector<string> binLabels;

};

#endif
