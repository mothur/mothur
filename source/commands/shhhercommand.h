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
#include "readcolumn.h"
#include "readmatrix.hpp"
#include "rabundvector.hpp"
#include "sabundvector.hpp"
#include "listvector.hpp"
#include "cluster.hpp"
#include "inputdata.h"
#include <cfloat>

//**********************************************************************************************************************

#define NUMBINS 1000
#define HOMOPS 10
#define MIN_COUNT 0.1
#define MIN_WEIGHT 0.1
#define MIN_TAU 0.0001
#define MIN_ITER 10
//**********************************************************************************************************************

class ShhherCommand : public Command {
	
public:
	ShhherCommand(string);
	~ShhherCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "shhh.flows";	}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Schloss PD, Gevers D, Westcott SL (2011).  Reducing the effects of PCR amplification and sequencing artifacts on 16S rRNA-based studies.  PLoS ONE.  6:e27310.\nQuince C, Lanzen A, Davenport RJ, Turnbaugh PJ (2011).  Removing noise from pyrosequenced amplicons.  BMC Bioinformatics  12:38.\nQuince C, Lanzén A, Curtis TP, Davenport RJ, Hall N, Head IM, Read LF, Sloan WT (2009).  Accurate determination of microbial diversity from 454 pyrosequencing data.  Nat. Methods 6:639.\nhttp://www.mothur.org/wiki/Shhh.flows"; }
	string getDescription()		{ return "shhh.flows"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }		
private:
	bool abort, large;
	string  flowFileName, flowFilesFileName, lookupFileName, compositeFASTAFileName, compositeNamesFileName;

	int maxIters, largeSize;
	float cutoff, sigma, minDelta;
	string flowOrder;
    
    vector<string> outputNames;
	vector<double> singleLookUp;
	vector<double> jointLookUp;
    vector<string> flowFileVector;
	
    vector<string> parseFlowFiles(string);
    int driver(vector<string>, string, string);
    int getFlowData(string, vector<string>&, vector<int>&, vector<short>&, map<string, int>&, int&);
    int getUniques(int, int, vector<short>&, vector<int>&, vector<int>&, vector<int>&, vector<int>&, vector<int>&, vector<double>&, vector<short>&);
    int flowDistParentFork(int, string, int, vector<int>&, vector<int>&, vector<int>&, vector<double>&, vector<short>&);
    float calcPairwiseDist(int, int, int, vector<int>&, vector<int>&, vector<double>&, vector<short>&);
    int createNamesFile(int, int, string, vector<string>&, vector<int>&, vector<int>&);
    int cluster(string, string, string);
    int getOTUData(int numSeqs, string,  vector<int>&, vector<int>&, vector<int>&, vector<vector<int> >&, vector<vector<int> >&, vector<int>&, vector<int>&,map<string, int>&);
    int calcCentroidsDriver(int numOTUs, vector<int>&, vector<int>&, vector<int>&, vector<short>&, vector<int>&, vector<double>&, vector<int>&, vector<short>&, vector<short>&, vector<int>&, int, vector<int>&);
    double getDistToCentroid(int, int, int, vector<short>&, vector<short>&, int);
    double getNewWeights(int, vector<int>&, vector<int>&, vector<double>&, vector<int>&, vector<double>&);
    
    double getLikelihood(int, int, vector<int>&, vector<int>&, vector<int>&, vector<int>&, vector<double>&, vector<double>&);
    int checkCentroids(int, vector<int>&, vector<double>&);
    void calcNewDistances(int, int, vector<int>& , vector<double>&,vector<double>& , vector<short>& change, vector<int>&,vector<vector<int> >&,	vector<double>&, vector<vector<int> >&, vector<int>&, vector<int>&, vector<short>&, vector<short>&, int, vector<int>&);
    int fill(int, vector<int>&, vector<int>&, vector<int>&, vector<int>&, vector<vector<int> >&, vector<vector<int> >&);
    void setOTUs(int, int, vector<int>&, vector<int>&, vector<int>&, vector<int>&,
                 vector<int>&, vector<double>&, vector<double>&, vector<vector<int> >&, vector<vector<int> >&);
    void writeQualities(int, int, string, vector<int>, vector<int>&, vector<int>&, vector<double>&, vector<short>&, vector<short>&, vector<int>&, vector<int>&, vector<string>&, vector<int>&, vector<vector<int> >&);
    void writeSequences(string, int, int, string, vector<int>, vector<short>&, vector<string>&, vector<vector<int> >&, vector<int>&);
    void writeNames(string, int, string, vector<int>, vector<string>&, vector<vector<int> >&, vector<int>&);
    void writeGroups(string, string, int, vector<string>&);
    void writeClusters(string, int, int, vector<int>, vector<int>&, vector<short>&, vector<string>&, vector<vector<int> >&, vector<int>&, vector<int>&, vector<short>&);
    
	void getSingleLookUp();
	void getJointLookUp();
    double getProbIntensity(int);
};

//**********************************************************************************************************************

#endif

