#ifndef SHAREDORDER_H
#define SHAREDORDER_H
/*
 *  sharedorder.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 12/9/08.
 *  Copyright 2008 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */
 
 /* This class is a child to datavector.  It represents OTU information at a certain distance. 
	It is similiar to an order vector except each member of data knows which group it belongs to.
	Each member of the internal container "data" represents is an individual which knows the OTU from which it came, 
	the group it is in and the abundance is equal to the OTU number.  */


#include "datavector.hpp"

struct individual {
		string group;
		int bin;
		int abundance;
		bool operator()(const individual& i1, const individual& i2) {
		return (i1.abundance > i2.abundance);
		}
};

struct individualFloat {
		string group;
		int bin;
		float abundance;
		bool operator()(const individual& i1, const individual& i2) {
		return (i1.abundance > i2.abundance);
		}
};


#include "sabundvector.hpp"
#include "rabundvector.hpp"
#include "sharedrabundvector.h"
#include "sharedsabundvector.h"
#include "globaldata.hpp"
#include "groupmap.h"
//#include "globaldata.hpp"

class GlobalData;

class SharedOrderVector : public DataVector {
	
public:
	SharedOrderVector();
//	SharedOrderVector(int ns, int nb=0, int mr=0)	: DataVector(), data(ns, -1), maxRank(0), numBins(0), numSeqs(0) {};
	SharedOrderVector(const SharedOrderVector& ov)	: DataVector(ov.label), data(ov.data), maxRank(ov.maxRank), numBins(ov.numBins), numSeqs(ov.numSeqs), needToUpdate(ov.needToUpdate) {if(needToUpdate == 1){	updateStats();}};

	SharedOrderVector(string, vector<individual>);
	SharedOrderVector(ifstream&);
	~SharedOrderVector(){};
	
	
	individual get(int);
	void resize(int);
	int size();
	void print(ostream&);
	vector<individual>::iterator begin();
	vector<individual>::iterator end();
	void push_back(int, int, string);  //OTU, abundance, group  MUST CALL UPDATE STATS AFTER PUSHBACK!!!
	void updateStats();

	int getNumBins();
	int getNumSeqs();
	int getMaxRank();
		
	RAbundVector getRAbundVector();
	SAbundVector getSAbundVector();
	OrderVector getOrderVector(map<string,int>*);
	SharedOrderVector getSharedOrderVector();
	SharedRAbundVector getSharedRAbundVector(string);  //get the sharedRabundvector for a sepecific group
	SharedSAbundVector getSharedSAbundVector(string);	//get the sharedSabundvector for a sepecific group
	vector<SharedRAbundVector*> getSharedRAbundVector(); //returns sharedRabundVectors for all the users groups
	
private:
	GlobalData* globaldata;
	GroupMap* groupmap;
	vector<individual>  data; 
	map< int, vector<individual> >::iterator it;
	int maxRank;
	int numBins;
	int numSeqs;
	bool needToUpdate;
	void set(int, int, int, string);	//index, OTU, abundance, group
	
};

#endif

