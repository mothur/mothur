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
	ShhherCommand();
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
	string outputDir, flowFileName, flowFilesFileName, lookupFileName, compositeFASTAFileName, compositeNamesFileName;

	int processors, maxIters, largeSize;
	float cutoff, sigma, minDelta;
	string flowOrder;
    
    vector<string> outputNames;
	vector<double> singleLookUp;
	vector<double> jointLookUp;
    vector<string> flowFileVector;
	
    vector<string> parseFlowFiles(string);
    int driver(vector<string>, string, string);
    int createProcesses(vector<string>);
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
	
	
#ifdef USE_MPI
	string flowDistMPI(int, int);
	void calcNewDistancesChildMPI(int, int, vector<int>&);

	int pid, ncpus;	
    
     void getFlowData();
     void getUniques();
     
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
    vector<int> nSeqsBreaks;
	vector<int> nOTUsBreaks;

#endif
	
};

/**************************************************************************************************
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct shhhFlowsData {
	int threadID, maxIters;
	float cutoff, sigma, minDelta;
	string flowOrder;
	vector<double> singleLookUp;
	vector<double> jointLookUp;
    vector<string> filenames;
    vector<string> outputNames;
    string thisCompositeFASTAFileName, thisCompositeNameFileName, outputDir;
    int start, stop;
	MothurOut* m;
	
	shhhFlowsData(){}
	shhhFlowsData(vector<string> f, string cf, string cn, string ou, string flor, vector<double> jl, vector<double> sl, MothurOut* mout, int st, int sp, float cut, float si, float mD, int mx, int tid) {
		filenames = f;
        thisCompositeFASTAFileName = cf;
        thisCompositeNameFileName = cn;
        outputDir = ou;
        flowOrder = flor;
		m = mout;
		start = st;
		stop = sp;
		cutoff= cut;
        sigma = si;
        minDelta = mD;
        maxIters = mx;
        jointLookUp = jl;
        singleLookUp = sl;
		threadID = tid;
	}
};

/**************************************************************************************************
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI ShhhFlowsThreadFunction(LPVOID lpParam){ 
	shhhFlowsData* pDataArray;
	pDataArray = (shhhFlowsData*)lpParam;
	
	try {
        
		for(int l=pDataArray->start;l<pDataArray->stop;l++){
			
			if (pDataArray->m->control_pressed) { break; }
			
			string flowFileName = pDataArray->filenames[l];
            
			pDataArray->m->mothurOut("\n>>>>>\tProcessing " + flowFileName + " (file " + toString(l+1) + " of " + toString(pDataArray->filenames.size()) + ")\t<<<<<\n");
			pDataArray->m->mothurOut("Reading flowgrams...\n");
			
            vector<string> seqNameVector;
            vector<int> lengths;
            vector<short> flowDataIntI;
            vector<double> flowDataPrI;
            map<string, int> nameMap;
            vector<short> uniqueFlowgrams;
            vector<int> uniqueCount;
            vector<int> mapSeqToUnique;
            vector<int> mapUniqueToSeq;
            vector<int> uniqueLengths;
            int numFlowCells;
            
            //int numSeqs = getFlowData(flowFileName, seqNameVector, lengths, flowDataIntI, nameMap, numFlowCells);
            /*****************************************************************************************************
            
            ifstream flowFile;
           // cout << "herethread " << flowFileName << '\t' << &flowFile << endl;
            pDataArray->m->openInputFile(flowFileName, flowFile);
            
           // cout << "herethread " << flowFileName << endl;
            string seqName;
            int currentNumFlowCells;
            float intensity;
                        
            flowFile >> numFlowCells;
            int index = 0;//pcluster
            while(!flowFile.eof()){
                
                if (pDataArray->m->control_pressed) { flowFile.close(); return 0; }
                
                flowFile >> seqName >> currentNumFlowCells;
                lengths.push_back(currentNumFlowCells);
             //  cout << "herethread " << seqName << endl;  
                seqNameVector.push_back(seqName);
                nameMap[seqName] = index++;//pcluster
                
                for(int i=0;i<numFlowCells;i++){
                    flowFile >> intensity;
                    if(intensity > 9.99)	{	intensity = 9.99;	}
                    int intI = int(100 * intensity + 0.0001);
                    flowDataIntI.push_back(intI);
                }
                pDataArray->m->gobble(flowFile);
            }
            flowFile.close();
            
            int numSeqs = seqNameVector.size();		
          //  cout << numSeqs << endl;   
            for(int i=0;i<numSeqs;i++){
                
                if (pDataArray->m->control_pressed) { return 0; }
                
                int iNumFlowCells = i * numFlowCells;
                for(int j=lengths[i];j<numFlowCells;j++){
                    flowDataIntI[iNumFlowCells + j] = 0;
                }
            }
          //  cout << "here" << endl; 
            /*****************************************************************************************************
	
			if (pDataArray->m->control_pressed) { return 0; }
			
			pDataArray->m->mothurOut("Identifying unique flowgrams...\n");
			//int numUniques = getUniques(numSeqs, numFlowCells, uniqueFlowgrams, uniqueCount, uniqueLengths, mapSeqToUnique, mapUniqueToSeq, lengths, flowDataPrI, flowDataIntI);
            /*****************************************************************************************************
            int numUniques = 0;
            uniqueFlowgrams.assign(numFlowCells * numSeqs, -1);
            uniqueCount.assign(numSeqs, 0);							//	anWeights
            uniqueLengths.assign(numSeqs, 0);
            mapSeqToUnique.assign(numSeqs, -1);
            mapUniqueToSeq.assign(numSeqs, -1);
            
            vector<short> uniqueFlowDataIntI(numFlowCells * numSeqs, -1);
            
            for(int i=0;i<numSeqs;i++){
                
                if (pDataArray->m->control_pressed) { return 0; }
                
                int index = 0;
                
                vector<short> current(numFlowCells);
                for(int j=0;j<numFlowCells;j++){
                    current[j] = short(((flowDataIntI[i * numFlowCells + j] + 50.0)/100.0));
                }
                
                for(int j=0;j<numUniques;j++){
                    int offset = j * numFlowCells;
                    bool toEnd = 1;
                    
                    int shorterLength;
                    if(lengths[i] < uniqueLengths[j])	{	shorterLength = lengths[i];			}
                    else								{	shorterLength = uniqueLengths[j];	}
                    
                    for(int k=0;k<shorterLength;k++){
                        if(current[k] != uniqueFlowgrams[offset + k]){
                            toEnd = 0;
                            break;
                        }
                    }
                    
                    if(toEnd){
                        mapSeqToUnique[i] = j;
                        uniqueCount[j]++;
                        index = j;
                        if(lengths[i] > uniqueLengths[j])	{	uniqueLengths[j] = lengths[i];	}
                        break;
                    }
                    index++;
                }
                
                if(index == numUniques){
                    uniqueLengths[numUniques] = lengths[i];
                    uniqueCount[numUniques] = 1;
                    mapSeqToUnique[i] = numUniques;//anMap
                    mapUniqueToSeq[numUniques] = i;//anF
                    
                    for(int k=0;k<numFlowCells;k++){
                        uniqueFlowgrams[numUniques * numFlowCells + k] = current[k];
                        uniqueFlowDataIntI[numUniques * numFlowCells + k] = flowDataIntI[i * numFlowCells + k];
                    }
                    
                    numUniques++;
                }
            }
            uniqueFlowDataIntI.resize(numFlowCells * numUniques);
            uniqueLengths.resize(numUniques);	
            
            flowDataPrI.resize(numSeqs * numFlowCells, 0);
            for(int i=0;i<flowDataPrI.size();i++)	{	
                if (pDataArray->m->control_pressed) { return 0; } 
                //flowDataPrI[i] = getProbIntensity(flowDataIntI[i]);	
                
                flowDataPrI[i] = 100000000; 
            
                for(int j=0;j<HOMOPS;j++){//loop signal strength
                    
                    if (pDataArray->m->control_pressed) { return 0; }
                    
                    float negLogProb = pDataArray->singleLookUp[j * NUMBINS + flowDataIntI[i]];
                    if(negLogProb < flowDataPrI[i])	{	flowDataPrI[i] = negLogProb; }
                }
            }            
            
            /*****************************************************************************************************
			
			if (pDataArray->m->control_pressed) { return 0; }
			
			pDataArray->m->mothurOut("Calculating distances between flowgrams...\n");
            string distFileName = flowFileName.substr(0,flowFileName.find_last_of('.')) + ".shhh.dist";
            unsigned long long begTime = time(NULL);
            double begClock = clock();
            
            //flowDistParentFork(numFlowCells, distFileName, numUniques, mapUniqueToSeq, mapSeqToUnique, lengths, flowDataPrI, flowDataIntI);	
            /*****************************************************************************************************
            ostringstream outStream;
            outStream.setf(ios::fixed, ios::floatfield);
            outStream.setf(ios::dec, ios::basefield);
            outStream.setf(ios::showpoint);
            outStream.precision(6);
            
            int thisbegTime = time(NULL);
            double thisbegClock = clock();
            
            for(int i=0;i<numUniques;i++){
                
                if (pDataArray->m->control_pressed) { break; }
                
                for(int j=0;j<i;j++){
                    //float flowDistance = calcPairwiseDist(numFlowCells, mapUniqueToSeq[i], mapUniqueToSeq[j], mapSeqToUnique, lengths, flowDataPrI, flowDataIntI);
                    /*****************************************************************************************************
                    int seqA = mapUniqueToSeq[i]; int seqB = mapUniqueToSeq[j];
                    int minLength = lengths[mapSeqToUnique[seqA]];
                    if(lengths[seqB] < minLength){	minLength = lengths[mapSeqToUnique[seqB]];	}
                    
                    int ANumFlowCells = seqA * numFlowCells;
                    int BNumFlowCells = seqB * numFlowCells;
                    
                    float flowDistance = 0;
                    
                    for(int k=0;k<minLength;k++){
                        
                        if (pDataArray->m->control_pressed) { break; }
                        
                        int flowAIntI = flowDataIntI[ANumFlowCells + k];
                        float flowAPrI = flowDataPrI[ANumFlowCells + k];
                        
                        int flowBIntI = flowDataIntI[BNumFlowCells + k];
                        float flowBPrI = flowDataPrI[BNumFlowCells + k];
                        flowDistance += pDataArray->jointLookUp[flowAIntI * NUMBINS + flowBIntI] - flowAPrI - flowBPrI;
                    }
                    
                    flowDistance /= (float) minLength;
                    /*****************************************************************************************************

                    if(flowDistance < 1e-6){
                        outStream << mapUniqueToSeq[i] << '\t' << mapUniqueToSeq[j] << '\t' << 0.000000 << endl;
                    }
                    else if(flowDistance <= pDataArray->cutoff){
                        outStream << mapUniqueToSeq[i] << '\t' << mapUniqueToSeq[j] << '\t' << flowDistance << endl;
                    }
                }
                if(i % 100 == 0){
                    pDataArray->m->mothurOut(toString(i) + "\t" + toString(time(NULL) - thisbegTime));
                    pDataArray->m->mothurOut("\t" + toString((clock()-thisbegClock)/CLOCKS_PER_SEC));
                    pDataArray->m->mothurOutEndLine();
                }
            }
            
            ofstream distFile(distFileName.c_str());
            distFile << outStream.str();		
            distFile.close();
            
            if (pDataArray->m->control_pressed) {}
            else {
                pDataArray->m->mothurOut(toString(numUniques-1) + "\t" + toString(time(NULL) - thisbegTime));
                pDataArray->m->mothurOut("\t" + toString((clock()-thisbegClock)/CLOCKS_PER_SEC));
                pDataArray->m->mothurOutEndLine();
            }
            /*****************************************************************************************************

            pDataArray->m->mothurOutEndLine();
            pDataArray->m->mothurOut("Total time: " + toString(time(NULL) - begTime) + '\t' + toString((clock() - begClock)/CLOCKS_PER_SEC) + '\n');
            
			string namesFileName = flowFileName.substr(0,flowFileName.find_last_of('.')) + ".shhh.names";
			//createNamesFile(numSeqs, numUniques, namesFileName, seqNameVector, mapSeqToUnique, mapUniqueToSeq);
            /*****************************************************************************************************
            vector<string> duplicateNames(numUniques, "");
            for(int i=0;i<numSeqs;i++){
                duplicateNames[mapSeqToUnique[i]] += seqNameVector[i] + ',';
            }
            
            ofstream nameFile;
            pDataArray->m->openOutputFile(namesFileName, nameFile);
            
            for(int i=0;i<numUniques;i++){
                if (pDataArray->m->control_pressed) { nameFile.close(); return 0; }
                nameFile << mapUniqueToSeq[i] << '\t' << duplicateNames[i].substr(0, duplicateNames[i].find_last_of(',')) << endl;
            }
            nameFile.close();
            /*****************************************************************************************************

			if (pDataArray->m->control_pressed) { return 0; }
			
			pDataArray->m->mothurOut("\nClustering flowgrams...\n");
            string listFileName = flowFileName.substr(0,flowFileName.find_last_of('.')) + ".shhh.list";
			//cluster(listFileName, distFileName, namesFileName);
            /*****************************************************************************************************
            ReadMatrix* read = new ReadColumnMatrix(distFileName); 	
            read->setCutoff(pDataArray->cutoff);
            
            NameAssignment* clusterNameMap = new NameAssignment(namesFileName);
            clusterNameMap->readMap();
            read->read(clusterNameMap);
            
            ListVector* list = read->getListVector();
            SparseMatrix* matrix = read->getMatrix();
            
            delete read; 
            delete clusterNameMap; 
            
            RAbundVector* rabund = new RAbundVector(list->getRAbundVector());
            
            Cluster* cluster = new CompleteLinkage(rabund, list, matrix, pDataArray->cutoff, "furthest"); 
            string tag = cluster->getTag();
            
            double clusterCutoff = pDataArray->cutoff;
            while (matrix->getSmallDist() <= clusterCutoff && matrix->getNNodes() > 0){
                
                if (pDataArray->m->control_pressed) { break; }
                
                cluster->update(clusterCutoff);
            }
            
            list->setLabel(toString(pDataArray->cutoff));
            
            ofstream listFileOut;
            pDataArray->m->openOutputFile(listFileName, listFileOut);
            list->print(listFileOut);
            listFileOut.close();
            
            delete matrix;	delete cluster;	delete rabund; delete list;
            /*****************************************************************************************************

			if (pDataArray->m->control_pressed) { return 0; }
            
            vector<int> otuData;
            vector<int> cumNumSeqs;
            vector<int> nSeqsPerOTU;
            vector<vector<int> > aaP;	//tMaster->aanP:	each row is a different otu / each col contains the sequence indices
            vector<vector<int> > aaI;	//tMaster->aanI:	that are in each otu - can't differentiate between aaP and aaI 
            vector<int> seqNumber;		//tMaster->anP:		the sequence id number sorted by OTU
            vector<int> seqIndex;		//tMaster->anI;		the index that corresponds to seqNumber
            
			
			//int numOTUs = getOTUData(numSeqs, listFileName, otuData, cumNumSeqs, nSeqsPerOTU, aaP, aaI, seqNumber, seqIndex, nameMap);
			/*****************************************************************************************************
            ifstream listFile;
            pDataArray->m->openInputFile(listFileName, listFile);
            string label;
            int numOTUs;
            
            listFile >> label >> numOTUs;
            
            otuData.assign(numSeqs, 0);
            cumNumSeqs.assign(numOTUs, 0);
            nSeqsPerOTU.assign(numOTUs, 0);
            aaP.clear();aaP.resize(numOTUs);
            
            seqNumber.clear();
            aaI.clear();
            seqIndex.clear();
            
            string singleOTU = "";
            
            for(int i=0;i<numOTUs;i++){
                
                if (pDataArray->m->control_pressed) { break; }
                
                listFile >> singleOTU;
                
                istringstream otuString(singleOTU);
                
                while(otuString){
                    
                    string seqName = "";
                    
                    for(int j=0;j<singleOTU.length();j++){
                        char letter = otuString.get();
                        
                        if(letter != ','){
                            seqName += letter;
                        }
                        else{
                            map<string,int>::iterator nmIt = nameMap.find(seqName);
                            int index = nmIt->second;
                            
                            nameMap.erase(nmIt);
                            
                            otuData[index] = i;
                            nSeqsPerOTU[i]++;
                            aaP[i].push_back(index);
                            seqName = "";
                        }
                    }
                    
                    map<string,int>::iterator nmIt = nameMap.find(seqName);
                    
                    int index = nmIt->second;
                    nameMap.erase(nmIt);
                    
                    otuData[index] = i;
                    nSeqsPerOTU[i]++;
                    aaP[i].push_back(index);	
                    
                    otuString.get();
                }
                
                sort(aaP[i].begin(), aaP[i].end());
                for(int j=0;j<nSeqsPerOTU[i];j++){
                    seqNumber.push_back(aaP[i][j]);
                }
                for(int j=nSeqsPerOTU[i];j<numSeqs;j++){
                    aaP[i].push_back(0);
                }
                
                
            }
            
            for(int i=1;i<numOTUs;i++){
                cumNumSeqs[i] = cumNumSeqs[i-1] + nSeqsPerOTU[i-1];
            }
            aaI = aaP;
            seqIndex = seqNumber;
            
            listFile.close();	    
            /*****************************************************************************************************

			if (pDataArray->m->control_pressed) { return 0; }
			
			pDataArray->m->mothurRemove(distFileName);
			pDataArray->m->mothurRemove(namesFileName);
			pDataArray->m->mothurRemove(listFileName);
			
            vector<double> dist;		//adDist - distance of sequences to centroids
            vector<short> change;		//did the centroid sequence change? 0 = no; 1 = yes
            vector<int> centroids;		//the representative flowgram for each cluster m
            vector<double> weight;
            vector<double> singleTau;	//tMaster->adTau:	1-D Tau vector (1xnumSeqs)
            vector<int> nSeqsBreaks;
            vector<int> nOTUsBreaks;
            
			dist.assign(numSeqs * numOTUs, 0);
            change.assign(numOTUs, 1);
            centroids.assign(numOTUs, -1);
            weight.assign(numOTUs, 0);
            singleTau.assign(numSeqs, 1.0);
            
            nSeqsBreaks.assign(2, 0);
            nOTUsBreaks.assign(2, 0);
            
            nSeqsBreaks[0] = 0;
            nSeqsBreaks[1] = numSeqs;
            nOTUsBreaks[1] = numOTUs;
			
			if (pDataArray->m->control_pressed) { break; }
			
			double maxDelta = 0;
			int iter = 0;
			
			begClock = clock();
			begTime = time(NULL);
            
			pDataArray->m->mothurOut("\nDenoising flowgrams...\n");
			pDataArray->m->mothurOut("iter\tmaxDelta\tnLL\t\tcycletime\n");
			
			while((pDataArray->maxIters == 0 && maxDelta > pDataArray->minDelta) || iter < MIN_ITER || (maxDelta > pDataArray->minDelta && iter < pDataArray->maxIters)){
				
				if (pDataArray->m->control_pressed) { break; }
				
				double cycClock = clock();
				unsigned long long cycTime = time(NULL);
				//fill(numOTUs, seqNumber, seqIndex, cumNumSeqs, nSeqsPerOTU, aaP, aaI);
                /*****************************************************************************************************
                int indexFill = 0;
                for(int i=0;i<numOTUs;i++){
                    
                    if (pDataArray->m->control_pressed) { return 0; }
                    
                    cumNumSeqs[i] = indexFill;
                    for(int j=0;j<nSeqsPerOTU[i];j++){
                        seqNumber[indexFill] = aaP[i][j];
                        seqIndex[indexFill] = aaI[i][j];
                        
                        indexFill++;
                    }
                }
                /*****************************************************************************************************

				
				if (pDataArray->m->control_pressed) { break; }
                
				//calcCentroidsDriver(numOTUs, cumNumSeqs, nSeqsPerOTU, seqIndex, change, centroids, singleTau, mapSeqToUnique, uniqueFlowgrams, flowDataIntI, lengths, numFlowCells, seqNumber);
                /*****************************************************************************************************
                for(int i=0;i<numOTUs;i++){
                    
                    if (pDataArray->m->control_pressed) { break; }
                    
                    double count = 0;
                    int position = 0;
                    int minFlowGram = 100000000;
                    double minFlowValue = 1e8;
                    change[i] = 0; //FALSE
                    
                    for(int j=0;j<nSeqsPerOTU[i];j++){
                        count += singleTau[seqNumber[cumNumSeqs[i] + j]];
                    }
                    
                    if(nSeqsPerOTU[i] > 0 && count > MIN_COUNT){
                        vector<double> adF(nSeqsPerOTU[i]);
                        vector<int> anL(nSeqsPerOTU[i]);
                        
                        for(int j=0;j<nSeqsPerOTU[i];j++){
                            int index = cumNumSeqs[i] + j;
                            int nI = seqIndex[index];
                            int nIU = mapSeqToUnique[nI];
                            
                            int k;
                            for(k=0;k<position;k++){
                                if(nIU == anL[k]){
                                    break;
                                }
                            }
                            if(k == position){
                                anL[position] = nIU;
                                adF[position] = 0.0000;
                                position++;
                            }						
                        }
                        
                        for(int j=0;j<nSeqsPerOTU[i];j++){
                            int index = cumNumSeqs[i] + j;
                            int nI = seqIndex[index];
                            
                            double tauValue = singleTau[seqNumber[index]];
                            
                            for(int k=0;k<position;k++){
                               // double dist = getDistToCentroid(anL[k], nI, lengths[nI], uniqueFlowgrams, flowDataIntI, numFlowCells);
                                /*****************************************************************************************************
                                int flowAValue = anL[k] * numFlowCells;
                                int flowBValue = nI * numFlowCells;
                                
                                double dist = 0;
                                
                                for(int l=0;l<lengths[nI];l++){
                                    dist += pDataArray->singleLookUp[uniqueFlowgrams[flowAValue] * NUMBINS + flowDataIntI[flowBValue]];
                                    flowAValue++;
                                    flowBValue++;
                                }
                                
                                dist = dist / (double)lengths[nI];
                                /*****************************************************************************************************
                                adF[k] += dist * tauValue;
                            }
                        }
                        
                        for(int j=0;j<position;j++){
                            if(adF[j] < minFlowValue){
                                minFlowGram = j;
                                minFlowValue = adF[j];
                            }
                        }
                        
                        if(centroids[i] != anL[minFlowGram]){
                            change[i] = 1;
                            centroids[i] = anL[minFlowGram];
                        }
                    }
                    else if(centroids[i] != -1){
                        change[i] = 1;
                        centroids[i] = -1;			
                    }
                }
                /*****************************************************************************************************

				if (pDataArray->m->control_pressed) { break; }
                
				//maxDelta = getNewWeights(numOTUs, cumNumSeqs, nSeqsPerOTU, singleTau, seqNumber, weight);  
                /*****************************************************************************************************
                double maxChange = 0;
                
                for(int i=0;i<numOTUs;i++){
                    
                    if (pDataArray->m->control_pressed) { break; }
                    
                    double difference = weight[i];
                    weight[i] = 0;
                    
                    for(int j=0;j<nSeqsPerOTU[i];j++){
                        int index = cumNumSeqs[i] + j;
                        double tauValue = singleTau[seqNumber[index]];
                        weight[i] += tauValue;
                    }
                    
                    difference = fabs(weight[i] - difference);
                    if(difference > maxChange){	maxChange = difference;	}
                }
                maxDelta = maxChange;
                /*****************************************************************************************************

                if (pDataArray->m->control_pressed) { break; }
                
				//double nLL = getLikelihood(numSeqs, numOTUs, nSeqsPerOTU, seqNumber, cumNumSeqs, seqIndex, dist, weight); 
                /*****************************************************************************************************
                vector<long double> P(numSeqs, 0);
                int effNumOTUs = 0;
                
                for(int i=0;i<numOTUs;i++){
                    if(weight[i] > MIN_WEIGHT){
                        effNumOTUs++;
                    }
                }
                
                string hold;
                for(int i=0;i<numOTUs;i++){
                    
                    if (pDataArray->m->control_pressed) { break; }
                    
                    for(int j=0;j<nSeqsPerOTU[i];j++){
                        int index = cumNumSeqs[i] + j;
                        int nI = seqIndex[index];
                        double singleDist = dist[seqNumber[index]];
                        
                        P[nI] += weight[i] * exp(-singleDist * pDataArray->sigma);
                    }
                }
                double nLL = 0.00;
                for(int i=0;i<numSeqs;i++){
                    if(P[i] == 0){	P[i] = DBL_EPSILON;	}
                    
                    nLL += -log(P[i]);
                }
                
                nLL = nLL -(double)numSeqs * log(pDataArray->sigma);
                /*****************************************************************************************************

                if (pDataArray->m->control_pressed) { break; }
                
				//checkCentroids(numOTUs, centroids, weight);
                /*****************************************************************************************************
                vector<int> unique(numOTUs, 1);
                
                for(int i=0;i<numOTUs;i++){
                    if(centroids[i] == -1 || weight[i] < MIN_WEIGHT){
                        unique[i] = -1;
                    }
                }
                
                for(int i=0;i<numOTUs;i++){
                    
                    if (pDataArray->m->control_pressed) { break; }
                    
                    if(unique[i] == 1){
                        for(int j=i+1;j<numOTUs;j++){
                            if(unique[j] == 1){
                                
                                if(centroids[j] == centroids[i]){
                                    unique[j] = 0;
                                    centroids[j] = -1;
                                    
                                    weight[i] += weight[j];
                                    weight[j] = 0.0;
                                }
                            }
                        }
                    }
                }
                /*****************************************************************************************************

				if (pDataArray->m->control_pressed) { break; }
				
				//calcNewDistances(numSeqs, numOTUs, nSeqsPerOTU,  dist, weight, change, centroids, aaP, singleTau, aaI, seqNumber, seqIndex, uniqueFlowgrams, flowDataIntI, numFlowCells, lengths);
                /*****************************************************************************************************
                int total = 0;
                vector<double> newTau(numOTUs,0);
                vector<double> norms(numSeqs, 0);
                nSeqsPerOTU.assign(numOTUs, 0);
                
                for(int i=0;i<numSeqs;i++){
                    
                    if (pDataArray->m->control_pressed) { break; }
                    
                    int indexOffset = i * numOTUs;
                    
                    double offset = 1e8;
                    
                    for(int j=0;j<numOTUs;j++){
                        
                        if(weight[j] > MIN_WEIGHT && change[j] == 1){
                            //dist[indexOffset + j] = getDistToCentroid(centroids[j], i, lengths[i], uniqueFlowgrams, flowDataIntI, numFlowCells);
                            /*****************************************************************************************************
                            int flowAValue = centroids[j] * numFlowCells;
                            int flowBValue = i * numFlowCells;
                            
                            double distTemp = 0;
                            
                            for(int l=0;l<lengths[i];l++){
                                distTemp += pDataArray->singleLookUp[uniqueFlowgrams[flowAValue] * NUMBINS + flowDataIntI[flowBValue]];
                                flowAValue++;
                                flowBValue++;
                            }
                            
                            dist[indexOffset + j] = distTemp / (double)lengths[i];
                            /*****************************************************************************************************

                        }
                        
                        if(weight[j] > MIN_WEIGHT && dist[indexOffset + j] < offset){
                            offset = dist[indexOffset + j];
                        }
                    }
                    
                    for(int j=0;j<numOTUs;j++){
                        if(weight[j] > MIN_WEIGHT){
                            newTau[j] = exp(pDataArray->sigma * (-dist[indexOffset + j] + offset)) * weight[j];
                            norms[i] += newTau[j];
                        }
                        else{
                            newTau[j] = 0.0;
                        }
                    }
                    
                    for(int j=0;j<numOTUs;j++){
                        newTau[j] /= norms[i];
                    }
                    
                    for(int j=0;j<numOTUs;j++){
                        if(newTau[j] > MIN_TAU){
                            
                            int oldTotal = total;
                            
                            total++;
                            
                            singleTau.resize(total, 0);
                            seqNumber.resize(total, 0);
                            seqIndex.resize(total, 0);
                            
                            singleTau[oldTotal] = newTau[j];
                            
                            aaP[j][nSeqsPerOTU[j]] = oldTotal;
                            aaI[j][nSeqsPerOTU[j]] = i;
                            nSeqsPerOTU[j]++;
                        }
                    }
                    
                }

                /*****************************************************************************************************

				if (pDataArray->m->control_pressed) { break; }
				
				iter++;
				
				pDataArray->m->mothurOut(toString(iter) + '\t' + toString(maxDelta) + '\t' + toString(nLL) + '\t' + toString(time(NULL) - cycTime) + '\t' + toString((clock() - cycClock)/(double)CLOCKS_PER_SEC) + '\n');
                
			}	
			
			if (pDataArray->m->control_pressed) { break; }
			
			pDataArray->m->mothurOut("\nFinalizing...\n");
			//fill(numOTUs, seqNumber, seqIndex, cumNumSeqs, nSeqsPerOTU, aaP, aaI);
            /*****************************************************************************************************
            int indexFill = 0;
            for(int i=0;i<numOTUs;i++){
                
                if (pDataArray->m->control_pressed) { return 0; }
                
                cumNumSeqs[i] = indexFill;
                for(int j=0;j<nSeqsPerOTU[i];j++){
                    seqNumber[indexFill] = aaP[i][j];
                    seqIndex[indexFill] = aaI[i][j];
                    
                    indexFill++;
                }
            }
            /*****************************************************************************************************

			if (pDataArray->m->control_pressed) { break; }
			
			//setOTUs(numOTUs, numSeqs, seqNumber, seqIndex, cumNumSeqs, nSeqsPerOTU, otuData, singleTau, dist, aaP, aaI);
            /*****************************************************************************************************
            vector<double> bigTauMatrix(numOTUs * numSeqs, 0.0000);
            
            for(int i=0;i<numOTUs;i++){
                
                if (pDataArray->m->control_pressed) { break; }
                
                for(int j=0;j<nSeqsPerOTU[i];j++){
                    int index = cumNumSeqs[i] + j;
                    double tauValue = singleTau[seqNumber[index]];
                    int sIndex = seqIndex[index];
                    bigTauMatrix[sIndex * numOTUs + i] = tauValue;				
                }
            }
            
            for(int i=0;i<numSeqs;i++){
                double maxTau = -1.0000;
                int maxOTU = -1;
                for(int j=0;j<numOTUs;j++){
                    if(bigTauMatrix[i * numOTUs + j] > maxTau){
                        maxTau = bigTauMatrix[i * numOTUs + j];
                        maxOTU = j;
                    }
                }
                
                otuData[i] = maxOTU;
            }
            
            nSeqsPerOTU.assign(numOTUs, 0);		
            
            for(int i=0;i<numSeqs;i++){
                int index = otuData[i];
                
                singleTau[i] = 1.0000;
                dist[i] = 0.0000;
                
                aaP[index][nSeqsPerOTU[index]] = i;
                aaI[index][nSeqsPerOTU[index]] = i;
                
                nSeqsPerOTU[index]++;
            }
            
            //fill(numOTUs, seqNumber, seqIndex, cumNumSeqs, nSeqsPerOTU, aaP, aaI);	
            /*****************************************************************************************************
            indexFill = 0;
            for(int i=0;i<numOTUs;i++){
                
                if (pDataArray->m->control_pressed) { return 0; }
                
                cumNumSeqs[i] = indexFill;
                for(int j=0;j<nSeqsPerOTU[i];j++){
                    seqNumber[indexFill] = aaP[i][j];
                    seqIndex[indexFill] = aaI[i][j];
                    
                    indexFill++;
                }
            }
            /*****************************************************************************************************/

            /*****************************************************************************************************

			if (pDataArray->m->control_pressed) { break; }
			
			vector<int> otuCounts(numOTUs, 0);
			for(int i=0;i<numSeqs;i++)	{	otuCounts[otuData[i]]++;	}
			
			//calcCentroidsDriver(numOTUs, cumNumSeqs, nSeqsPerOTU, seqIndex, change, centroids, singleTau, mapSeqToUnique, uniqueFlowgrams, flowDataIntI, lengths, numFlowCells, seqNumber);	
            /*****************************************************************************************************
            for(int i=0;i<numOTUs;i++){
                
                if (pDataArray->m->control_pressed) { break; }
                
                double count = 0;
                int position = 0;
                int minFlowGram = 100000000;
                double minFlowValue = 1e8;
                change[i] = 0; //FALSE
                
                for(int j=0;j<nSeqsPerOTU[i];j++){
                    count += singleTau[seqNumber[cumNumSeqs[i] + j]];
                }
                
                if(nSeqsPerOTU[i] > 0 && count > MIN_COUNT){
                    vector<double> adF(nSeqsPerOTU[i]);
                    vector<int> anL(nSeqsPerOTU[i]);
                    
                    for(int j=0;j<nSeqsPerOTU[i];j++){
                        int index = cumNumSeqs[i] + j;
                        int nI = seqIndex[index];
                        int nIU = mapSeqToUnique[nI];
                        
                        int k;
                        for(k=0;k<position;k++){
                            if(nIU == anL[k]){
                                break;
                            }
                        }
                        if(k == position){
                            anL[position] = nIU;
                            adF[position] = 0.0000;
                            position++;
                        }						
                    }
                    
                    for(int j=0;j<nSeqsPerOTU[i];j++){
                        int index = cumNumSeqs[i] + j;
                        int nI = seqIndex[index];
                        
                        double tauValue = singleTau[seqNumber[index]];
                        
                        for(int k=0;k<position;k++){
                            // double dist = getDistToCentroid(anL[k], nI, lengths[nI], uniqueFlowgrams, flowDataIntI, numFlowCells);
                            /*****************************************************************************************************
                            int flowAValue = anL[k] * numFlowCells;
                            int flowBValue = nI * numFlowCells;
                            
                            double dist = 0;
                            
                            for(int l=0;l<lengths[nI];l++){
                                dist += pDataArray->singleLookUp[uniqueFlowgrams[flowAValue] * NUMBINS + flowDataIntI[flowBValue]];
                                flowAValue++;
                                flowBValue++;
                            }
                            
                            dist = dist / (double)lengths[nI];
                            /*****************************************************************************************************
                            adF[k] += dist * tauValue;
                        }
                    }
                    
                    for(int j=0;j<position;j++){
                        if(adF[j] < minFlowValue){
                            minFlowGram = j;
                            minFlowValue = adF[j];
                        }
                    }
                    
                    if(centroids[i] != anL[minFlowGram]){
                        change[i] = 1;
                        centroids[i] = anL[minFlowGram];
                    }
                }
                else if(centroids[i] != -1){
                    change[i] = 1;
                    centroids[i] = -1;			
                }
            }

            /*****************************************************************************************************

            if (pDataArray->m->control_pressed) { break; }
            
			//writeQualities(numOTUs, numFlowCells, flowFileName, otuCounts, nSeqsPerOTU, seqNumber, singleTau, flowDataIntI, uniqueFlowgrams, cumNumSeqs, mapUniqueToSeq, seqNameVector, centroids, aaI); 
            if (pDataArray->m->control_pressed) { break; }
            /*****************************************************************************************************
            string thisOutputDir = pDataArray->outputDir;
            if (pDataArray->outputDir == "") {  thisOutputDir += pDataArray->m->hasPath(flowFileName);  }
            string qualityFileName = thisOutputDir + pDataArray->m->getRootName(pDataArray->m->getSimpleName(flowFileName)) + "shhh.qual";
            
            ofstream qualityFile;
            pDataArray->m->openOutputFile(qualityFileName, qualityFile);
            
            qualityFile.setf(ios::fixed, ios::floatfield);
            qualityFile.setf(ios::showpoint);
            qualityFile << setprecision(6);
            
            vector<vector<int> > qualities(numOTUs);
            vector<double> pr(HOMOPS, 0);
            
            
            for(int i=0;i<numOTUs;i++){
                
                if (pDataArray->m->control_pressed) { break; }
                
                int index = 0;
                int base = 0;
                
                if(nSeqsPerOTU[i] > 0){
                    qualities[i].assign(1024, -1);
                    
                    while(index < numFlowCells){
                        double maxPrValue = 1e8;
                        short maxPrIndex = -1;
                        double count = 0.0000;
                        
                        pr.assign(HOMOPS, 0);
                        
                        for(int j=0;j<nSeqsPerOTU[i];j++){
                            int lIndex = cumNumSeqs[i] + j;
                            double tauValue = singleTau[seqNumber[lIndex]];
                            int sequenceIndex = aaI[i][j];
                            short intensity = flowDataIntI[sequenceIndex * numFlowCells + index];
                            
                            count += tauValue;
                            
                            for(int s=0;s<HOMOPS;s++){
                                pr[s] += tauValue * pDataArray->singleLookUp[s * NUMBINS + intensity];
                            }
                        }
                        
                        maxPrIndex = uniqueFlowgrams[centroids[i] * numFlowCells + index];
                        maxPrValue = pr[maxPrIndex];
                        
                        if(count > MIN_COUNT){
                            double U = 0.0000;
                            double norm = 0.0000;
                            
                            for(int s=0;s<HOMOPS;s++){
                                norm += exp(-(pr[s] - maxPrValue));
                            }
                            
                            for(int s=1;s<=maxPrIndex;s++){
                                int value = 0;
                                double temp = 0.0000;
                                
                                U += exp(-(pr[s-1]-maxPrValue))/norm;
                                
                                if(U>0.00){
                                    temp = log10(U);
                                }
                                else{
                                    temp = -10.1;
                                }
                                temp = floor(-10 * temp);
                                value = (int)floor(temp);
                                if(value > 100){	value = 100;	}
                                
                                qualities[i][base] = (int)value;
                                base++;
                            }
                        }
                        
                        index++;
                    }
                }
                
                
                if(otuCounts[i] > 0){
                    qualityFile << '>' << seqNameVector[mapUniqueToSeq[i]] << endl;
                    
                    int j=4;	//need to get past the first four bases
                    while(qualities[i][j] != -1){
                        qualityFile << qualities[i][j] << ' ';
                        j++;
                    }
                    qualityFile << endl;
                }
            }
            qualityFile.close();
            pDataArray->outputNames.push_back(qualityFileName);
            /*****************************************************************************************************

           // writeSequences(thisCompositeFASTAFileName, numOTUs, numFlowCells, flowFileName, otuCounts, uniqueFlowgrams, seqNameVector, aaI, centroids);
            if (pDataArray->m->control_pressed) { break; }
            /*****************************************************************************************************
            thisOutputDir = pDataArray->outputDir;
            if (pDataArray->outputDir == "") {  thisOutputDir += pDataArray->m->hasPath(flowFileName);  }
            string fastaFileName = thisOutputDir + pDataArray->m->getRootName(pDataArray->m->getSimpleName(flowFileName)) + "shhh.fasta";
            ofstream fastaFile;
            pDataArray->m->openOutputFile(fastaFileName, fastaFile);
            
            vector<string> names(numOTUs, "");
            
            for(int i=0;i<numOTUs;i++){
                
                if (pDataArray->m->control_pressed) { break; }
                
                int index = centroids[i];
                
                if(otuCounts[i] > 0){
                    fastaFile << '>' << seqNameVector[aaI[i][0]] << endl;
                    
                    string newSeq = "";
                    
                    for(int j=0;j<numFlowCells;j++){
                        
                        char base = pDataArray->flowOrder[j % 4];
                        for(int k=0;k<uniqueFlowgrams[index * numFlowCells + j];k++){
                            newSeq += base;
                        }
                    }
                    
                    fastaFile << newSeq.substr(4) << endl;
                }
            }
            fastaFile.close();
            
            pDataArray->outputNames.push_back(fastaFileName);
            
            if(pDataArray->thisCompositeFASTAFileName != ""){
                pDataArray->m->appendFiles(fastaFileName, pDataArray->thisCompositeFASTAFileName);
            }

            /*****************************************************************************************************

            //writeNames(thisCompositeNamesFileName, numOTUs, flowFileName, otuCounts, seqNameVector, aaI, nSeqsPerOTU);				
            if (pDataArray->m->control_pressed) { break; }
            /*****************************************************************************************************
            thisOutputDir = pDataArray->outputDir;
            if (pDataArray->outputDir == "") {  thisOutputDir += pDataArray->m->hasPath(flowFileName);  }
            string nameFileName = thisOutputDir + pDataArray->m->getRootName(pDataArray->m->getSimpleName(flowFileName)) + "shhh.names";
            ofstream nameFileOut;
            pDataArray->m->openOutputFile(nameFileName, nameFileOut);
            
            for(int i=0;i<numOTUs;i++){
                
                if (pDataArray->m->control_pressed) { break; }
                
                if(otuCounts[i] > 0){
                    nameFileOut << seqNameVector[aaI[i][0]] << '\t' << seqNameVector[aaI[i][0]];
                    
                    for(int j=1;j<nSeqsPerOTU[i];j++){
                        nameFileOut << ',' << seqNameVector[aaI[i][j]];
                    }
                    
                    nameFileOut << endl;
                }
            }
            nameFileOut.close();
            pDataArray->outputNames.push_back(nameFileName);
            
            
            if(pDataArray->thisCompositeNameFileName != ""){
                pDataArray->m->appendFiles(nameFileName, pDataArray->thisCompositeNameFileName);
            }		
            /*****************************************************************************************************

            //writeClusters(flowFileName, numOTUs, numFlowCells,otuCounts, centroids, uniqueFlowgrams, seqNameVector, aaI, nSeqsPerOTU, lengths, flowDataIntI);			
            if (pDataArray->m->control_pressed) { break; }
            /*****************************************************************************************************
            thisOutputDir = pDataArray->outputDir;
            if (pDataArray->outputDir == "") {  thisOutputDir += pDataArray->m->hasPath(flowFileName);  }
            string otuCountsFileName = thisOutputDir + pDataArray->m->getRootName(pDataArray->m->getSimpleName(flowFileName)) + "shhh.counts";
            ofstream otuCountsFile;
            pDataArray->m->openOutputFile(otuCountsFileName, otuCountsFile);
            
            string bases = pDataArray->flowOrder;
            
            for(int i=0;i<numOTUs;i++){
                
                if (pDataArray->m->control_pressed) {
                    break;
                }
                //output the translated version of the centroid sequence for the otu
                if(otuCounts[i] > 0){
                    int index = centroids[i];
                    
                    otuCountsFile << "ideal\t";
                    for(int j=8;j<numFlowCells;j++){
                        char base = bases[j % 4];
                        for(int s=0;s<uniqueFlowgrams[index * numFlowCells + j];s++){
                            otuCountsFile << base;
                        }
                    }
                    otuCountsFile << endl;
                    
                    for(int j=0;j<nSeqsPerOTU[i];j++){
                        int sequence = aaI[i][j];
                        otuCountsFile << seqNameVector[sequence] << '\t';
                        
                        string newSeq = "";
                        
                        for(int k=0;k<lengths[sequence];k++){
                            char base = bases[k % 4];
                            int freq = int(0.01 * (double)flowDataIntI[sequence * numFlowCells + k] + 0.5);
                            
                            for(int s=0;s<freq;s++){
                                newSeq += base;
                                //otuCountsFile << base;
                            }
                        }
                        otuCountsFile << newSeq.substr(4) << endl;
                    }
                    otuCountsFile << endl;
                }
            }
            otuCountsFile.close();
            pDataArray->outputNames.push_back(otuCountsFileName)
            /*****************************************************************************************************

            //writeGroups(flowFileName, numSeqs, seqNameVector);						
            if (pDataArray->m->control_pressed) { break; }
            /*****************************************************************************************************
            thisOutputDir = pDataArray->outputDir;
            if (pDataArray->outputDir == "") {  thisOutputDir += pDataArray->m->hasPath(flowFileName);  }
            string fileRoot = thisOutputDir + pDataArray->m->getRootName(pDataArray->m->getSimpleName(flowFileName));
            string groupFileName = fileRoot + "shhh.groups";
            ofstream groupFile;
            pDataArray->m->openOutputFile(groupFileName, groupFile);
            
            for(int i=0;i<numSeqs;i++){
                if (pDataArray->m->control_pressed) { break; }
                groupFile << seqNameVector[i] << '\t' << fileRoot << endl;
            }
            groupFile.close();
            pDataArray->outputNames.push_back(groupFileName);
            /*****************************************************************************************************

            pDataArray->m->mothurOut("Total time to process " + flowFileName + ":\t" + toString(time(NULL) - begTime) + '\t' + toString((clock() - begClock)/(double)CLOCKS_PER_SEC) + '\n');
		}
		
        if (pDataArray->m->control_pressed) { for (int i = 0; i < pDataArray->outputNames.size(); i++) { pDataArray->m->mothurRemove(pDataArray->outputNames[i]); } return 0; }
        
        return 0;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "ShhherCommand", "ShhhFlowsThreadFunction");
		exit(1);
	}
} 
#endif
*/

#endif

