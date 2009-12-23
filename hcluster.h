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
#include "nameassignment.hpp"

class RAbundVector;
class ListVector;

/***********************************************************************/
class HCluster {
	
public:
	HCluster(RAbundVector*, ListVector*, string);
	~HCluster(){};
    void update(int, int, float);
	void setMapWanted(bool m); 
	map<string, int> getSeqtoBin()  {  return seq2Bin;	}
	vector<seqDist> getSeqs(ifstream&, NameAssignment*, float);

protected:	
	void clusterBins();
	void clusterNames();
	int getUpmostParent(int);
	int makeActive();
	void printInfo();
	void updateArrayandLinkTable();
	void updateMap();
		
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
	map<string, int> seq2Bin;
	bool mapWanted, exitedBreak;
	seqDist next;
	string method;
	
};

/***********************************************************************/







#endif


