#ifndef SHAREDRABUNDFLOATVECTOR_H
#define SHAREDRABUNDFLOATVECTOR_H

/*
 *  sharedrabundfloatvector.h
 *  Mothur
 *
 *  Created by westcott on 8/18/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


#include "datavector.hpp"
#include "sharedordervector.h"
#include "sharedsabundvector.h"
#include "rabundvector.hpp"
//#include "groupmap.h"

/* This class is a child to datavector.  It represents OTU information at a certain distance. 
	It is similiar to an rabundvector except each member of data knows which group it belongs to.
	Each member of the internal container "data" is a struct of type individualFloat. 
	An individual which knows the OTU from which it came, 
	the group it is in and its abundance.  */


class SharedRAbundFloatVector : public DataVector {
	
public:
	SharedRAbundFloatVector();
	SharedRAbundFloatVector(int);
	SharedRAbundFloatVector(const SharedRAbundFloatVector& bv) : DataVector(bv), data(bv.data), maxRank(bv.maxRank), numBins(bv.numBins), numSeqs(bv.numSeqs), group(bv.group), index(bv.index){};
    SharedRAbundFloatVector(ifstream&);
	~SharedRAbundFloatVector();

	int getNumBins();		
	float getNumSeqs();							
	float getMaxRank();
	string getGroup();
	void setGroup(string);
	int getGroupIndex();
	void setGroupIndex(int);

	void set(int, float, string);			//OTU, abundance, groupname
	individualFloat get(int);
	vector <individual> getData();
	float getAbundance(int);
    vector<float> getAbundances();
	void push_front(float, int, string); //abundance, otu, groupname
	void insert(float, int, string); //abundance, otu, groupname
	void push_back(float, string);  //abundance, groupname
	void pop_back();
	void resize(int);
	void clear();
	int size();
	
	void print(ostream&);
	void printHeaders(ostream&);
		
	RAbundVector getRAbundVector();
	SAbundVector getSAbundVector();
	OrderVector getOrderVector(map<string,int>*);
	//SharedOrderVector getSharedOrderVector();
	//SharedSAbundVector getSharedSAbundVector();
	//SharedRAbundVector getSharedRAbundVector();
	SharedRAbundFloatVector getSharedRAbundFloatVector();
	vector<SharedRAbundFloatVector*> getSharedRAbundFloatVectors();
	
private:
    int eliminateZeroOTUS();
	vector<individualFloat>  data; 
	vector<SharedRAbundFloatVector*> lookup;
	//GlobalData* globaldata;
	//GroupMap* groupmap;
	float maxRank;
	int numBins;
	float numSeqs;
	string group;
	int index;	
	
	
};


#endif


