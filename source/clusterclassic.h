#ifndef CLUSTERCLASSIC_H
#define CLUSTERCLASSIC_H


#include "mothurout.h"
#include "listvector.hpp"
#include "rabundvector.hpp"
#include "nameassignment.hpp"
#include "counttable.h"

/*
 *  clusterclassic.h
 *  Mothur
 *
 *  Created by westcott on 10/29/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */


class ClusterClassic {
	
public:
	ClusterClassic(float, string, bool);
	int readPhylipFile(string, NameAssignment*);
    int readPhylipFile(string, CountTable*);
	void update(double&);
	double getSmallDist() { return smallDist; }	
	int getNSeqs() { return nseqs; }	
	ListVector* getListVector() { return list; }
	RAbundVector* getRAbundVector() { return rabund; }		
	string getTag() { return tag; }
	void setMapWanted(bool m);  
	map<string, int> getSeqtoBin()  {  return seq2Bin;	}

private:	
	double getSmallCell();
	void clusterBins();
	void clusterNames();
	void updateMap();
	void print();
	
	struct colDist {
		int col;
		int row;
		float dist;
		colDist(int r, int c, double d) : row(r), col(c), dist(d) {}
	};
	
	RAbundVector* rabund;
	ListVector* list;
	vector< vector<float> > dMatrix;	
    Utils util;
	
	int smallRow;
	int smallCol, nseqs;
	double smallDist;
	bool mapWanted, sim;
	double cutoff, aboveCutoff;
	map<string, int> seq2Bin;
	string method, tag;
	
	MothurOut* m;
};

#endif


