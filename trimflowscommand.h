#ifndef TRIMFLOWSCOMMAND_H
#define TRIMFLOWSCOMMAND_H

/*
 *  trimflowscommand.h
 *  Mothur
 *
 *  Created by Pat Schloss on 12/22/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "sequence.hpp"
#include "flowdata.h"
#include "groupmap.h"
#include "trimoligos.h"

class TrimFlowsCommand : public Command {
public:
	TrimFlowsCommand(string);
	TrimFlowsCommand();
	~TrimFlowsCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "trim.flows";	}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Trim.flows"; }
	string getDescription()		{ return "trim.flows"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	bool abort;

	struct linePair {
		unsigned long long start;
		unsigned long long end;
		linePair(unsigned long long i, unsigned long long j) : start(i), end(j) {}
	};
	int comboStarts;
	vector<int> processIDS;   //processid
	vector<linePair*> lines;

	vector<unsigned long long> getFlowFileBreaks();
	int createProcessesCreateTrim(string, string, string, string, vector<vector<string> >); 
	int driverCreateTrim(string, string, string, string, vector<vector<string> >, linePair*);
    string reverseOligo(string);
    
	vector<string> outputNames;
	set<string> filesToRemove;
	
	void getOligos(vector<vector<string> >&);		//a rewrite of what is in trimseqscommand.h
	
	bool allFiles;
	int processors;
	int numFPrimers, numRPrimers;
	int maxFlows, minFlows, minLength, maxLength, maxHomoP, tdiffs, bdiffs, pdiffs, sdiffs, ldiffs, numLinkers, numSpacers;
	int numFlows;
	float signal, noise;
	bool fasta;
	string flowOrder;	
	
	string flowFileName, oligoFileName, outputDir;

	map<string, int> barcodes;
	map<string, int> primers;
	vector<string> revPrimer;
    vector<string> linker;
    vector<string> spacer;

	vector<string> primerNameVector;	//needed here?
	vector<string> barcodeNameVector;	//needed here?

	map<string, int> combos;			//needed here?
	map<string, int> groupToIndex;		//needed here?
	
};

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct trimFlowData {
	string flowFileName; 
	string trimFlowFileName; 
	string scrapFlowFileName;
	string fastaFileName;
	string flowOrder;
	vector<vector<string> > barcodePrimerComboFileNames;
	map<string, int> barcodes;
	map<string, int> primers;
	vector<string> revPrimer;
	bool fasta, allFiles;
	unsigned long long start;
	unsigned long long end;
	MothurOut* m;
	float signal, noise;
	int numFlows, maxFlows, minFlows, maxHomoP, tdiffs, bdiffs, pdiffs, threadID, count;
	
	trimFlowData(){}
	trimFlowData(string ff, string tf, string sf, string f, string fo, vector<vector<string> > bfn, map<string, int> bar, map<string, int> pri, vector<string> rev, bool fa, bool al, unsigned long long st, unsigned long long en, MothurOut* mout, float sig, float n, int numF, int maxF, int minF, int maxH, int td, int bd, int pd, int tid) {
		flowFileName = ff;
		trimFlowFileName = tf;
		scrapFlowFileName = sf;
		fastaFileName = f;
		flowOrder = fo;
		barcodePrimerComboFileNames = bfn;
		barcodes = bar;
		primers = pri;
		revPrimer = rev;
		fasta = fa;
		allFiles = al;
		start = st;
		end = en;
		m = mout;
		signal = sig;
		noise = n;
		numFlows = numF;
		maxFlows = maxF;
		minFlows = minF;
		maxHomoP = maxH;
		tdiffs = td;
		bdiffs = bd;
		pdiffs = pd;
		threadID = tid;
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyTrimFlowThreadFunction(LPVOID lpParam){ 
	trimFlowData* pDataArray;
	pDataArray = (trimFlowData*)lpParam;
	
	try {
		ofstream trimFlowFile;
		pDataArray->m->openOutputFile(pDataArray->trimFlowFileName, trimFlowFile);
		trimFlowFile.setf(ios::fixed, ios::floatfield); trimFlowFile.setf(ios::showpoint);
		
		ofstream scrapFlowFile;
		pDataArray->m->openOutputFile(pDataArray->scrapFlowFileName, scrapFlowFile);
		scrapFlowFile.setf(ios::fixed, ios::floatfield); scrapFlowFile.setf(ios::showpoint);
		
		ofstream fastaFile;
		if(pDataArray->fasta){	pDataArray->m->openOutputFile(pDataArray->fastaFileName, fastaFile);	}
		
		ifstream flowFile;
		pDataArray->m->openInputFile(pDataArray->flowFileName, flowFile);
		
		flowFile.seekg(pDataArray->start);
		
		if(pDataArray->start == 0){
			flowFile >> pDataArray->numFlows; pDataArray->m->gobble(flowFile);
			scrapFlowFile << pDataArray->maxFlows << endl;
			trimFlowFile << pDataArray->maxFlows << endl;
			if(pDataArray->allFiles){
				for(int i=0;i<pDataArray->barcodePrimerComboFileNames.size();i++){
					for(int j=0;j<pDataArray->barcodePrimerComboFileNames[0].size();j++){
						ofstream temp;
						pDataArray->m->openOutputFile(pDataArray->barcodePrimerComboFileNames[i][j], temp);
						temp << pDataArray->maxFlows << endl;
						temp.close();
					}
				}			
			}
		}
		
		FlowData flowData(pDataArray->numFlows, pDataArray->signal, pDataArray->noise, pDataArray->maxHomoP, pDataArray->flowOrder);
		cout << " thread flowdata address " <<  &flowData  << '\t' << &flowFile << endl;
		TrimOligos trimOligos(pDataArray->pdiffs, pDataArray->bdiffs, pDataArray->primers, pDataArray->barcodes, pDataArray->revPrimer);
		
		pDataArray->count = pDataArray->end;
		cout << pDataArray->threadID << '\t' << pDataArray->count << endl;
		int count = 0;
		for(int i = 0; i < pDataArray->end; i++){ //end is the number of sequences to process
			
			if (pDataArray->m->control_pressed) {  break; }
			cout << pDataArray->threadID << '\t' << count << endl;
			int success = 1;
			int currentSeqDiffs = 0;
			string trashCode = "";
				
			flowData.getNext(flowFile);
			cout << "thread good bit " << flowFile.good() << endl;
			flowData.capFlows(pDataArray->maxFlows);	
			
			Sequence currSeq = flowData.getSequence();
			if(!flowData.hasMinFlows(pDataArray->minFlows)){	//screen to see if sequence is of a minimum number of flows
				success = 0;
				trashCode += 'l';
			}
			
			int primerIndex = 0;
			int barcodeIndex = 0;
			
			if(pDataArray->barcodes.size() != 0){
				success = trimOligos.stripBarcode(currSeq, barcodeIndex);
				if(success > pDataArray->bdiffs)		{	trashCode += 'b';	}
				else{ currentSeqDiffs += success;  }
			}
			
			if(pDataArray->primers.size() != 0){
				success = trimOligos.stripForward(currSeq, primerIndex);
				if(success > pDataArray->pdiffs)		{	trashCode += 'f';	}
				else{ currentSeqDiffs += success;  }
			}
			
			if (currentSeqDiffs > pDataArray->tdiffs)	{	trashCode += 't';   }
			
			if(pDataArray->revPrimer.size() != 0){
				success = trimOligos.stripReverse(currSeq);
				if(!success)				{	trashCode += 'r';	}
			}
			
			if(trashCode.length() == 0){
				
				flowData.printFlows(trimFlowFile);
				
				if(pDataArray->fasta)	{	currSeq.printSequence(fastaFile);	}
				
				if(pDataArray->allFiles){
					ofstream output;
					pDataArray->m->openOutputFileAppend(pDataArray->barcodePrimerComboFileNames[barcodeIndex][primerIndex], output);
					output.setf(ios::fixed, ios::floatfield); trimFlowFile.setf(ios::showpoint);
					
					flowData.printFlows(output);
					output.close();
				}				
			}
			else{
				flowData.printFlows(scrapFlowFile, trashCode);
			}
			
			count++;
				cout << pDataArray->threadID << '\t' << currSeq.getName() << endl;		
			//report progress
			if((count) % 10000 == 0){	pDataArray->m->mothurOut(toString(count)); pDataArray->m->mothurOutEndLine();		}
			
		}
		//report progress
		if((count) % 10000 != 0){	pDataArray->m->mothurOut(toString(count)); pDataArray->m->mothurOutEndLine();		}
		
		trimFlowFile.close();
		scrapFlowFile.close();
		flowFile.close();
		if(pDataArray->fasta){	fastaFile.close();	}
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "TrimFlowsCommand", "MyTrimFlowsThreadFunction");
		exit(1);
	}
} 
#endif


#endif
