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
	HCluster(RAbundVector*, ListVector*, string, string, NameAssignment*, float);
	~HCluster(){};
    bool update(int, int, float);
	void setMapWanted(bool m); 
	map<string, int> getSeqtoBin()  {  return seq2Bin;	}
	vector<seqDist> getSeqs();

protected:	
	void clusterBins();
	void clusterNames();
	int getUpmostParent(int);
	int makeActive();
	void printInfo();
	void updateArrayandLinkTable();
	void updateMap();
	vector<seqDist> getSeqsFNNN();
	vector<seqDist> getSeqsAN();
	void combineFile();
	void processFile();
	//seqDist getNextDist(char*, int&, int);
		
	RAbundVector* rabund;
	ListVector* list;
	NameAssignment* nameMap;
	
	vector<clusterNode> clusterArray;
	
	//note: the nearest and average neighbor method do not use the link table or active links
	vector< map<int, int> > linkTable;  // vector of maps - linkTable[1][6] = 2  would mean sequence in spot 1 has 2 links with sequence in 6
	map<int, int> activeLinks;  //maps sequence to index in linkTable
	map<int, int>::iterator it;
	map<int, int>::iterator itActive;
	map<int, int>::iterator it2Active;
	map<int, int>::iterator it2;
	
	int numSeqs;
	int smallRow;
	int smallCol;
	float smallDist, cutoff;
	map<string, int> seq2Bin;
	bool mapWanted, exitedBreak;
	seqDist next;
	string method, distfile;
	ifstream filehandle;
	
	vector<seqDist> mergedMin;
	string partialDist;
	
	
};

/***********************************************************************/







#endif


