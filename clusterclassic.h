#ifndef CLUSTER_H
#define CLUSTER_H


#include "mothur.h"
#include "mothurout.h"
#include "listvector.hpp"
#include "rabundvector.hpp"
#include "nameassignment.hpp"

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
	void update(double&);
	double getSmallDist() { return smallDist; }	
	int getNSeqs() { return nseqs; }	
	ListVector* getListVector() { return list; }
	RAbundVector* getRAbundVector() { return rabund; }		
	string getTag();
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
		double dist;
		colDist(int r, int c, double d) : row(r), col(c), dist(d) {}
	};
	
	RAbundVector* rabund;
	ListVector* list;
	vector< vector<double> > dMatrix;	
	//vector<colDist> rowSmallDists;
	
	int smallRow;
	int smallCol, nseqs;
	double smallDist;
	bool mapWanted, sim;
	double cutoff, aboveCutoff;
	map<string, int> seq2Bin;
	string method;
	
	MothurOut* m;
};

#endif


