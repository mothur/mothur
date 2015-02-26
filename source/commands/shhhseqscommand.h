#ifndef SHHHSEQSCOMMAND_H
#define SHHHSEQSCOMMAND_H

/*
 *  shhhseqscommand.h
 *  Mothur
 *
 *  Created by westcott on 11/8/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "command.hpp"
#include "myseqdist.h"
#include "seqnoise.h"
#include "sequenceparser.h"
#include "deconvolutecommand.h"
#include "clustercommand.h"

//**********************************************************************************************************************

class ShhhSeqsCommand : public Command {
	
public:
	ShhhSeqsCommand(string);
	ShhhSeqsCommand();
	~ShhhSeqsCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "shhh.seqs";	}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Schloss PD, Gevers D, Westcott SL (2011).  Reducing the effects of PCR amplification and sequencing artifacts on 16S rRNA-based studies.  PLoS ONE.  6:e27310.\nQuince C, Lanzen A, Davenport RJ, Turnbaugh PJ (2011).  Removing noise from pyrosequenced amplicons.  BMC Bioinformatics  12:38.\nhttp://www.mothur.org/wiki/Shhh.seqs"; }
	string getDescription()		{ return "shhh.seqs"; }
	
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }		
	
private:
	
	struct linePair {
		int start;
		int end;
		linePair(int i, int j) : start(i), end(j) {}
	};
	
	bool abort;
	string outputDir, fastafile, namefile, groupfile;
	int processors;
	double sigma;
	vector<string> outputNames;
	
	int readData(correctDist*, seqNoise&, vector<string>&, vector<string>&, vector<string>&, vector<int>&);
	int loadData(correctDist*, seqNoise&, vector<string>&, vector<string>&, vector<string>&, vector<int>&, map<string, string>&, vector<Sequence>&);

	int driver(seqNoise&, vector<string>&, vector<string>&, vector<string>&, vector<int>&, string, string, string, string);
	vector<string> driverGroups(SequenceParser&, string, string, string, int, int, vector<string>);
	vector<string> createProcessesGroups(SequenceParser&, string, string, string, vector<string>);
	int deconvoluteResults(string, string);



};

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct shhhseqsData {
	string fastafile; 
	string namefile; 
	string groupfile;
	string newFName, newNName, newMName;
	MothurOut* m;
	int start;
	int end;
	int sigma, threadID, count;
	vector<string> groups;
	vector<string> mapfileNames;
	
	shhhseqsData(){}
	shhhseqsData(string f, string n, string g, string nff,  string nnf, string nmf, vector<string> gr, MothurOut* mout, int st, int en, int s, int tid) {
		fastafile = f;
		namefile = n;
		groupfile = g;
		newFName = nff;
		newNName = nnf;
		newMName = nmf;
		m = mout;
		start = st;
		end = en;
		sigma = s;
		threadID = tid;
		groups = gr;
        count=0;
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyShhhSeqsThreadFunction(LPVOID lpParam){ 
	shhhseqsData* pDataArray;
	pDataArray = (shhhseqsData*)lpParam;
	
	try {
		
		//parse fasta and name file by group
		SequenceParser parser(pDataArray->groupfile, pDataArray->fastafile, pDataArray->namefile);
						
		//precluster each group
		for (int k = pDataArray->start; k < pDataArray->end; k++) {
			
            pDataArray->count++;
            
			int start = time(NULL);
			
			if (pDataArray->m->control_pressed) {  return 0; }
			
			pDataArray->m->mothurOutEndLine(); pDataArray->m->mothurOut("Processing group " + pDataArray->groups[k] + ":"); pDataArray->m->mothurOutEndLine();

			map<string, string> thisNameMap;
			thisNameMap = parser.getNameMap(pDataArray->groups[k]); 
			vector<Sequence> thisSeqs = parser.getSeqs(pDataArray->groups[k]);
			
			if (pDataArray->m->control_pressed) {  return 0; }
			
			vector<string> sequences;
			vector<string> uniqueNames;
			vector<string> redundantNames;
			vector<int> seqFreq;
			
			seqNoise noise;
			correctDist* correct = new correctDist(1); //we use one processor since we already split up the work load.
			
			//load this groups info in order
			//loadData(correct, noise, sequences, uniqueNames, redundantNames, seqFreq, thisNameMap, thisSeqs);
			//////////////////////////////////////////////////////////////////////////////////////////////////
			bool error = false;
			map<string, string>::iterator it;
			
			for (int i = 0; i < thisSeqs.size(); i++) {
				
				if (pDataArray->m->control_pressed) { return 0; }
				
				if (thisSeqs[i].getName() != "") {
					correct->addSeq(thisSeqs[i].getName(), thisSeqs[i].getAligned());
					
					it = thisNameMap.find(thisSeqs[i].getName());
					if (it != thisNameMap.end()) {
						noise.addSeq(thisSeqs[i].getAligned(), sequences);
						noise.addRedundantName(it->first, it->second, uniqueNames, redundantNames, seqFreq);
					}else {
						pDataArray->m->mothurOut("[ERROR]: " + thisSeqs[i].getName() + " is in your fasta file and not in your namefile, please correct.");
						error = true;
					}
				}
			}
			
			if (error) { return 0; }
			//////////////////////////////////////////////////////////////////////////////////////////////////
			
			if (pDataArray->m->control_pressed) { return 0; }
			
			//calc distances for cluster
			string distFileName = pDataArray->m->getRootName(pDataArray->m->getSimpleName(pDataArray->fastafile)) + pDataArray->groups[k] + ".shhh.dist";
			correct->execute(distFileName);
			delete correct;
			
			if (pDataArray->m->control_pressed) { pDataArray->m->mothurRemove(distFileName); return 0; }
			
			//driver(noise, sequences, uniqueNames, redundantNames, seqFreq, distFileName, newFFile+groups[i], newNFile+groups[i], newMFile+groups[i]+".map"); 
			///////////////////////////////////////////////////////////////////////////////////////////////////
			double cutOff = 0.08;
			int minIter = 10;
			int maxIter = 1000;
			double minDelta = 1e-6;
			int numIters = 0;
			double maxDelta = 1e6;
			int numSeqs = sequences.size();
			
			//run cluster command
			string inputString = "phylip=" + distFileName + ", method=furthest, cutoff=0.08";
			pDataArray->m->mothurOut("/******************************************/"); pDataArray->m->mothurOutEndLine(); 
			pDataArray->m->mothurOut("Running command: cluster(" + inputString + ")"); pDataArray->m->mothurOutEndLine(); 
			
			Command* clusterCommand = new ClusterCommand(inputString);
			clusterCommand->execute();
			
			map<string, vector<string> > filenames = clusterCommand->getOutputFiles();
			string listFileName = filenames["list"][0];
			string rabundFileName = filenames["rabund"][0]; pDataArray->m->mothurRemove(rabundFileName);
			string sabundFileName = filenames["sabund"][0]; pDataArray->m->mothurRemove(sabundFileName);
			
			delete clusterCommand;
			pDataArray->m->mothurOut("/******************************************/"); pDataArray->m->mothurOutEndLine(); 
			
			if (pDataArray->m->control_pressed) { pDataArray->m->mothurRemove(distFileName); pDataArray->m->mothurRemove(listFileName); return 0; }
			
			vector<double> distances(numSeqs * numSeqs);
			noise.getDistanceData(distFileName, distances);
			pDataArray->m->mothurRemove(distFileName); 
			if (pDataArray->m->control_pressed) { pDataArray->m->mothurRemove(listFileName); return 0; }
			
			vector<int> otuData(numSeqs);
			vector<int> otuFreq;
			vector<vector<int> > otuBySeqLookUp;
			noise.getListData(listFileName, cutOff, otuData, otuFreq, otuBySeqLookUp);
			pDataArray->m->mothurRemove(listFileName);
			if (pDataArray->m->control_pressed) { return 0; }
			
			int numOTUs = otuFreq.size();
			
			vector<double> weights(numOTUs, 0);
			vector<int> change(numOTUs, 1);
			vector<int> centroids(numOTUs, -1);
			vector<int> cumCount(numOTUs, 0);
			
			vector<double> tau(numSeqs, 1);
			vector<int> anP(numSeqs, 0);
			vector<int> anI(numSeqs, 0);
			vector<int> anN(numSeqs, 0);
			vector<vector<int> > aanI = otuBySeqLookUp;
			
			while(numIters < minIter || ((maxDelta > minDelta) && (numIters < maxIter))){
				
				if (pDataArray->m->control_pressed) { return 0; }
				
				noise.updateOTUCountData(otuFreq, otuBySeqLookUp, aanI, anP, anI, cumCount); if (pDataArray->m->control_pressed) { return 0; }
				maxDelta = noise.calcNewWeights(weights, seqFreq, anI, cumCount, anP, otuFreq, tau);  if (pDataArray->m->control_pressed) { return 0; }
				
				noise.calcCentroids(anI, anP, change, centroids, cumCount, distances, seqFreq, otuFreq, tau); if (pDataArray->m->control_pressed) { return 0; }
				noise.checkCentroids(weights, centroids); if (pDataArray->m->control_pressed) { return 0; }
				
				otuFreq.assign(numOTUs, 0);
				
				int total = 0;
				
				for(int i=0;i<numSeqs;i++){
					if (pDataArray->m->control_pressed) { return 0; }
					
					double offset = 1e6;
					double norm = 0.0000;
					double minWeight = 0.1;
					vector<double> currentTau(numOTUs);
					
					for(int j=0;j<numOTUs;j++){
						if (pDataArray->m->control_pressed) { return 0; }
						if(weights[j] > minWeight && distances[i * numSeqs+centroids[j]] < offset){
							offset = distances[i * numSeqs+centroids[j]];
						}
					}
					
					for(int j=0;j<numOTUs;j++){
						if (pDataArray->m->control_pressed) { return 0; }
						if(weights[j] > minWeight){
							currentTau[j] = exp(pDataArray->sigma * (-distances[(i * numSeqs + centroids[j])] + offset)) * weights[j];
							norm += currentTau[j];
						}
						else{
							currentTau[j] = 0.0000;
						}
					}			
					
					for(int j=0;j<numOTUs;j++){
						if (pDataArray->m->control_pressed) { return 0; }
						currentTau[j] /= norm;
					}
					
					for(int j=0;j<numOTUs;j++){
						if (pDataArray->m->control_pressed) { return 0; }
						
						if(currentTau[j] > 1.0e-4){
							int oldTotal = total;
							total++;
							
							tau.resize(oldTotal+1);
							tau[oldTotal] = currentTau[j];
							otuBySeqLookUp[j][otuFreq[j]] = oldTotal;
							aanI[j][otuFreq[j]] = i;
							otuFreq[j]++;
							
						}
					}
					
					anP.resize(total);
					anI.resize(total);
				}
				
				numIters++;
			}
			
			noise.updateOTUCountData(otuFreq, otuBySeqLookUp, aanI, anP, anI, cumCount);  if (pDataArray->m->control_pressed) { return 0; }
			
			vector<double> percentage(numSeqs);
			noise.setUpOTUData(otuData, percentage, cumCount, tau, otuFreq, anP, anI);  if (pDataArray->m->control_pressed) { return 0; }
			noise.finishOTUData(otuData, otuFreq, anP, anI, cumCount, otuBySeqLookUp, aanI, tau);  if (pDataArray->m->control_pressed) { return 0; }
			
			change.assign(numOTUs, 1);
			noise.calcCentroids(anI, anP, change, centroids, cumCount, distances, seqFreq, otuFreq, tau); if (pDataArray->m->control_pressed) { return 0; }
			
			
			vector<int> finalTau(numOTUs, 0);
			for(int i=0;i<numSeqs;i++){
				if (pDataArray->m->control_pressed) { return 0; }
				finalTau[otuData[i]] += int(seqFreq[i]);
			}
			
			noise.writeOutput(pDataArray->newFName+pDataArray->groups[k], pDataArray->newNName+pDataArray->groups[k], pDataArray->newMName+pDataArray->groups[k]+".map", finalTau, centroids, otuData, sequences, uniqueNames, redundantNames, seqFreq, distances);
			
			///////////////////////////////////////////////////////////////////////////////////////////////////
			
			if (pDataArray->m->control_pressed) { return 0; }
			
			pDataArray->m->appendFiles(pDataArray->newFName+pDataArray->groups[k], pDataArray->newFName); pDataArray->m->mothurRemove(pDataArray->newFName+pDataArray->groups[k]);
			pDataArray->m->appendFiles(pDataArray->newNName+pDataArray->groups[k], pDataArray->newNName); pDataArray->m->mothurRemove(pDataArray->newNName+pDataArray->groups[k]);
			pDataArray->mapfileNames.push_back(pDataArray->newMName+pDataArray->groups[k]+".map");
			
			pDataArray->m->mothurOut("It took " + toString(time(NULL) - start) + " secs to process group " + pDataArray->groups[k] + "."); pDataArray->m->mothurOutEndLine(); 
		}
		
		return 0;
		
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "ShhhSeqsCommand", "MyShhhSeqsThreadFunction");
		exit(1);
	}
} 
#endif
/**************************************************************************************************/

#endif
