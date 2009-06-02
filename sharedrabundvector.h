#ifndef SHAREDRABUNDVECTOR_H
#define SHAREDRABUNDVECTOR_H

/*
 *  sharedrabundvector.h
 *  Dotur
 *
 *  Created by Sarah Westcott on 12/5/08.
 *  Copyright 2008 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "datavector.hpp"
#include "sharedordervector.h"
#include "sharedsabundvector.h"
#include "rabundvector.hpp"
#include "groupmap.h"

/* This class is a child to datavector.  It represents OTU information at a certain distance. 
	It is similiar to an rabundvector except each member of data knows which group it belongs to.
	Each member of the internal container "data" is a struct of type individual. 
	An individual which knows the OTU from which it came, 
	the group it is in and its abundance.  */

class GlobalData;

class SharedRAbundVector : public DataVector {
	
public:
	SharedRAbundVector();
	SharedRAbundVector(int);
	//SharedRAbundVector(string, vector<int>);
	SharedRAbundVector(const SharedRAbundVector& bv) : DataVector(bv), data(bv.data), maxRank(bv.maxRank), numBins(bv.numBins), numSeqs(bv.numSeqs){};
    SharedRAbundVector(ifstream&);
	~SharedRAbundVector(){};

	int getNumBins();		
	int getNumSeqs();							
	int getMaxRank();
	string getGroup();
	void setGroup(string);
	int getGroupIndex();
	void setGroupIndex(int);								

	void set(int, int, string);			//OTU, abundance, groupname
	void setData(vector <individual>);
	individual get(int);
	vector <individual> getData();
	int getAbundance(int);
	int numNZ();
	void sortD();  //Sorts the data in descending order.
	void push_front(int, int, string); //abundance, otu, groupname
	void insert(int, int, string); //abundance, otu, groupname
	void push_back(int, int, string);  //abundance, otu, groupname
	void pop_back();
	void resize(int);
	int size();
	vector<individual>::reverse_iterator rbegin();
	vector<individual>::reverse_iterator rend();
	
	void print(ostream&);
		
	RAbundVector getRAbundVector();
	RAbundVector getRAbundVector2();
	SAbundVector getSAbundVector();
	OrderVector getOrderVector(map<string,int>*);
	SharedOrderVector getSharedOrderVector();
	SharedSAbundVector getSharedSAbundVector();
	SharedRAbundVector getSharedRAbundVector();
	vector<SharedRAbundVector*> getSharedRAbundVectors();
	
private:
	vector<individual>  data; 
	vector<SharedRAbundVector*> lookup;
	GlobalData* globaldata;
	GroupMap* groupmap;
	int maxRank;
	int numBins;
	int numSeqs;
	string group;
	int index;	
};


#endif

