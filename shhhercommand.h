#ifndef SHHHER_H
#define SHHHER_H

/*
 *  shhher.h
 *  Mothur
 *
 *  Created by Pat Schloss on 12/27/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "globaldata.hpp"

class ShhherCommand : public Command {
	
public:
	ShhherCommand(string);
	ShhherCommand();
	~ShhherCommand();
	vector<string> getRequiredParameters();
	vector<string> getValidParameters();
	vector<string> getRequiredFiles();
	map<string, vector<string> > getOutputFiles() { return outputTypes; }
	int execute();
	void help();
	
private:
	GlobalData* globaldata;
	
	int abort;
	map<string, vector<string> > outputTypes;
	
	string outputDir, flowFileName, flowFilesFileName, lookupFileName, compositeFASTAFileName;

	int processors, maxIters;
	float cutoff, sigma, minDelta;
	string flowOrder;
	
	vector<int> nSeqsBreaks;
	vector<int> nOTUsBreaks;
	vector<double> singleLookUp;
	vector<double> jointLookUp;
	
	vector<string> seqNameVector;
	vector<int> lengths;
	vector<short> flowDataIntI;
	vector<double> flowDataPrI;
	map<string, int> nameMap;
	vector<int> otuData;
	vector<int> cumNumSeqs;
	vector<int> nSeqsPerOTU;
	vector<vector<int> > aaP;	//tMaster->aanP:	each row is a different otu / each col contains the sequence indices
	vector<vector<int> > aaI;	//tMaster->aanI:	that are in each otu - can't differentiate between aaP and aaI 
	vector<int> seqNumber;		//tMaster->anP:		the sequence id number sorted by OTU
	vector<int> seqIndex;		//tMaster->anI;		the index that corresponds to seqNumber
	vector<double> dist;		//adDist - distance of sequences to centroids
	vector<short> change;		//did the centroid sequence change? 0 = no; 1 = yes
	vector<int> centroids;		//the representative flowgram for each cluster m
	vector<double> weight;
	vector<double> singleTau;	//tMaster->adTau:	1-D Tau vector (1xnumSeqs)
	vector<short> uniqueFlowgrams;
	vector<int> uniqueCount;
	vector<int> mapSeqToUnique;
	vector<int> mapUniqueToSeq;
	vector<int> uniqueLengths;
	
	int numSeqs, numUniques, numOTUs, numFlowCells;
	
	void getSingleLookUp();
	void getJointLookUp();
	void getFlowData();
	void getUniques();
	double getProbIntensity(int);
	float calcPairwiseDist(int, int);
	void flowDistParentFork(string, int, int);
	
	string createDistFile(int);
	string createNamesFile();
	string cluster(string, string);
	
	void getOTUData(string);
	void initPyroCluster();
	void fill();
	void calcCentroids();
	void calcCentroidsDriver(int, int);
	double getDistToCentroid(int, int, int);
	double getNewWeights();
	double getLikelihood();
	void checkCentroids();
	void calcNewDistances();
	void calcNewDistancesParent(int, int);
	void calcNewDistancesChild(int, int, vector<int>&, vector<int>&, vector<double>&);


	void setOTUs();
	void writeQualities(vector<int>);
	void writeSequences(vector<int>);
	void writeNames(vector<int>);
	void writeGroups();
	void writeClusters(vector<int>);

	
#ifdef USE_MPI
	string flowDistMPI(int, int);
	void calcNewDistancesChildMPI(int, int, vector<int>&);

	int pid, ncpus;	
#endif
	
};


#endif

