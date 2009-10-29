#ifndef HCLUSTER_H
#define HCLUSTER_H

/*
 *  hcluster.h
 *  Mothur
 *
 *  Created by westcott on 10/13/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */


#include "mothur.h"

class RAbundVector;
class ListVector;
/************************************************************/
struct clusterNode {
	int numSeq;
	int parent;
	int smallChild; //used to make linkTable work with list and rabund
	clusterNode(int num, int par, int kid) : numSeq(num), parent(par), smallChild(kid) {};
};

/***********************************************************************/
class HCluster {
	
public:
	HCluster(RAbundVector*, ListVector*);
    bool update(int, int, float);
	//string getTag();

protected:	
	void clusterBins();
	void clusterNames();
	int getUpmostParent(int);
	int makeActive();
	void printInfo();
	void updateArrayandLinkTable();
		
	RAbundVector* rabund;
	ListVector* list;
	
	vector<clusterNode> clusterArray;
	vector< map<int, int> > linkTable;  // vector of maps - linkTable[1][6] = 2  would mean sequence in spot 1 has 2 links with sequence in 6
	map<int, int> activeLinks;  //maps sequence to index in linkTable
	map<int, int>::iterator it;
	map<int, int>::iterator it2;
	
	int numSeqs;
	
	int smallRow;
	int smallCol;
	float smallDist;
	
};

/***********************************************************************/







#endif


