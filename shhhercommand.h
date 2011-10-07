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
	ShhherCommand();
	~ShhherCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "shhh.flows";	}
	string getCommandCategory()		{ return "Sequence Processing";		}
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Shhh.flows"; }
	string getDescription()		{ return "shhh.flows"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }		
private:
	
	int abort;
	
	string outputDir, flowFileName, flowFilesFileName, lookupFileName, compositeFASTAFileName, compositeNamesFileName;

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

	vector<string> outputNames;

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

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
typedef struct flowDistParentForkData {
	string distFileName; 
	vector<int> mapUniqueToSeq;
	vector<int> mapSeqToUnique;
	vector<int> lengths;
	vector<short> flowDataIntI;
	vector<double> flowDataPrI;
	vector<double> jointLookUp;
	MothurOut* m;
	int threadID, startSeq, stopSeq, numFlowCells;
	float cutoff;
	
	flowDistParentForkData(){}
	flowDistParentForkData(string d, vector<int> mapU, vector<int> mapS, vector<int> l, vector<short> flowD, vector<double> flowDa, vector<double> j, MothurOut* mout, int st, int sp, int n, float cut, int tid) {
		distFileName = d;
		mapUniqueToSeq = mapU;
		mapSeqToUnique = mapS;
		lengths = l;
		flowDataIntI = flowD;
		flowDataPrI = flowDa;
		jointLookUp = j;
		m = mout;
		startSeq = st;
		stopSeq = sp;
		numFlowCells = n;
		cutoff= cut;
		threadID = tid;
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
#else
static DWORD WINAPI MyflowDistParentForkThreadFunction(LPVOID lpParam){ 
	flowDistParentForkData* pDataArray;
	pDataArray = (flowDistParentForkData*)lpParam;
	
	try {
		ostringstream outStream;
		outStream.setf(ios::fixed, ios::floatfield);
		outStream.setf(ios::dec, ios::basefield);
		outStream.setf(ios::showpoint);
		outStream.precision(6);
		
		int begTime = time(NULL);
		double begClock = clock();
		string tempOut = "start and end = " + toString(pDataArray->startSeq) +'\t' + toString(pDataArray->stopSeq) + "-";
		cout << tempOut << endl;
	
		for(int i=pDataArray->startSeq;i<pDataArray->stopSeq;i++){
			
			if (pDataArray->m->control_pressed) { break; }
			cout << "thread i = " << i << endl;
			for(int j=0;j<i;j++){
				
				cout << "thread j = " << j << endl;
				float flowDistance = 0.0;
				////////////////// calcPairwiseDist ///////////////////
				//needed because this is a static global function that can't see the classes internal functions
				
				int minLength = pDataArray->lengths[pDataArray->mapSeqToUnique[pDataArray->mapUniqueToSeq[i]]];
				if(pDataArray->lengths[pDataArray->mapUniqueToSeq[j]] < minLength){	minLength = pDataArray->lengths[pDataArray->mapSeqToUnique[pDataArray->mapUniqueToSeq[j]]];	}
				
				int ANumFlowCells = pDataArray->mapUniqueToSeq[i] * pDataArray->numFlowCells;
				int BNumFlowCells = pDataArray->mapUniqueToSeq[j] * pDataArray->numFlowCells;
				
				for(int k=0;k<minLength;k++){
					
					if (pDataArray->m->control_pressed) { break; }
					
					int flowAIntI = pDataArray->flowDataIntI[ANumFlowCells + k];
					float flowAPrI = pDataArray->flowDataPrI[ANumFlowCells + k];
					
					int flowBIntI = pDataArray->flowDataIntI[BNumFlowCells + k];
					float flowBPrI = pDataArray->flowDataPrI[BNumFlowCells + k];
					flowDistance += pDataArray->jointLookUp[flowAIntI * NUMBINS + flowBIntI] - flowAPrI - flowBPrI;
				}
				
				flowDistance /= (float) minLength;
				//cout << flowDistance << endl;
				////////////////// end of calcPairwiseDist ///////////////////
								
				if(flowDistance < 1e-6){
					outStream << pDataArray->mapUniqueToSeq[i] << '\t' << pDataArray->mapUniqueToSeq[j] << '\t' << 0.000000 << endl;
				}
				else if(flowDistance <= pDataArray->cutoff){
					outStream << pDataArray->mapUniqueToSeq[i] << '\t' << pDataArray->mapUniqueToSeq[j] << '\t' << flowDistance << endl;
				}
			}
			if(i % 100 == 0){
				pDataArray->m->mothurOut(toString(i) + "\t" + toString(time(NULL) - begTime));
				pDataArray->m->mothurOut("\t" + toString((clock()-begClock)/CLOCKS_PER_SEC));
				pDataArray->m->mothurOutEndLine();
			}
		}
		
		ofstream distFile(pDataArray->distFileName.c_str());
		distFile << outStream.str();		
		distFile.close();
		
		if (pDataArray->m->control_pressed) {}
		else {
			pDataArray->m->mothurOut(toString(pDataArray->stopSeq-1) + "\t" + toString(time(NULL) - begTime));
			pDataArray->m->mothurOut("\t" + toString((clock()-begClock)/CLOCKS_PER_SEC));
			pDataArray->m->mothurOutEndLine();
		}		
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "ShhherCommand", "MyflowDistParentForkThreadFunction");
		exit(1);
	}
} 
#endif


#endif

