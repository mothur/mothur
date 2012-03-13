#ifndef SCREENSEQSCOMMAND_H
#define SCREENSEQSCOMMAND_H

/*
 *  screenseqscommand.h
 *  Mothur
 *
 *  Created by Pat Schloss on 6/3/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */
#include "mothur.h"
#include "command.hpp"
#include "sequence.hpp"

class ScreenSeqsCommand : public Command {
	
public:
	ScreenSeqsCommand(string);
	ScreenSeqsCommand();
	~ScreenSeqsCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "screen.seqs";				}
	string getCommandCategory()		{ return "Sequence Processing";		}
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Screen.seqs"; }
	string getDescription()		{ return "enables you to keep sequences that fulfill certain user defined criteria"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:

	struct linePair {
		unsigned long long start;
		unsigned long long end;
		linePair(unsigned long long i, unsigned long long j) : start(i), end(j) {}
	};

	vector<linePair> lines;

	int screenNameGroupFile(set<string>);
	int screenGroupFile(set<string>);
	int screenAlignReport(set<string>);
	int screenQual(set<string>);
	int screenTaxonomy(set<string>);
	
	int driver(linePair, string, string, string, set<string>&);
	int createProcesses(string, string, string, set<string>&);
	
	#ifdef USE_MPI
	int driverMPI(int, int, MPI_File&, MPI_File&, MPI_File&, vector<unsigned long long>&, set<string>&);
	#endif

	bool abort;
	string fastafile, namefile, groupfile, alignreport, outputDir, qualfile, taxonomy;
	int startPos, endPos, maxAmbig, maxHomoP, minLength, maxLength, processors, criteria;
	vector<string> outputNames;
	vector<string> optimize;
	map<string, int> nameMap;
	int readNames();
	
	int getSummary(vector<unsigned long long>&);
	int createProcessesCreateSummary(vector<int>&, vector<int>&, vector<int>&, vector<int>&, vector<int>&, string);
	int driverCreateSummary(vector<int>&, vector<int>&, vector<int>&, vector<int>&, vector<int>&, string, linePair);	
};

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct sumData {
	vector<int> startPosition;
	vector<int> endPosition;
	vector<int> seqLength; 
	vector<int> ambigBases; 
	vector<int> longHomoPolymer; 
	string filename, namefile; 
	unsigned long long start;
	unsigned long long end;
	int count;
	MothurOut* m;
	map<string, int> nameMap;
	
	
	sumData(){}
	sumData(string f, MothurOut* mout, unsigned long long st, unsigned long long en, string nf, map<string, int> nam) {
		filename = f;
        namefile = nf;
		m = mout;
		start = st;
		end = en;
		nameMap = nam;
		count = 0;
	}
};
/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct sumScreenData {
    int startPos, endPos, maxAmbig, maxHomoP, minLength, maxLength;
	unsigned long long start;
	unsigned long long end;
	int count;
	MothurOut* m;
	string goodFName, badAccnosFName, filename;
    set<string>* badSeqNames;
	
	
	sumScreenData(){}
	sumScreenData(int s, int e, int a, int h, int minl, int maxl, string f, MothurOut* mout, unsigned long long st, unsigned long long en, string gf, string bf, set<string>* bn) {
		startPos = s;
		endPos = e;
		minLength = minl;
        maxLength = maxl;
		maxAmbig = a;
		maxHomoP = h;
		filename = f;
        goodFName = gf;
        badAccnosFName = bf;
		m = mout;
		start = st;
		end = en;
		badSeqNames = bn;
		count = 0;
	}
};


/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MySumThreadFunction(LPVOID lpParam){ 
	sumData* pDataArray;
	pDataArray = (sumData*)lpParam;
	
	try {
		ifstream in;
		pDataArray->m->openInputFile(pDataArray->filename, in);
        
		//print header if you are process 0
		if ((pDataArray->start == 0) || (pDataArray->start == 1)) {
			in.seekg(0);
		}else { //this accounts for the difference in line endings. 
			in.seekg(pDataArray->start-1); pDataArray->m->gobble(in); 
		}
		
		pDataArray->count = pDataArray->end;
		for(int i = 0; i < pDataArray->end; i++){ //end is the number of sequences to process
			
			if (pDataArray->m->control_pressed) { in.close();  pDataArray->count = 1; return 1; }
			
			Sequence current(in); pDataArray->m->gobble(in); 
			
			if (current.getName() != "") {
				
				int num = 1;
				if (pDataArray->namefile != "") {
					//make sure this sequence is in the namefile, else error 
					map<string, int>::iterator it = pDataArray->nameMap.find(current.getName());
					
					if (it == pDataArray->nameMap.end()) { pDataArray->m->mothurOut("[ERROR]: " + current.getName() + " is not in your namefile, please correct."); pDataArray->m->mothurOutEndLine(); pDataArray->m->control_pressed = true; }
					else { num = it->second; }
				}
				
				//for each sequence this sequence represents
				for (int i = 0; i < num; i++) {
					pDataArray->startPosition.push_back(current.getStartPos());
					pDataArray->endPosition.push_back(current.getEndPos());
					pDataArray->seqLength.push_back(current.getNumBases());
					pDataArray->ambigBases.push_back(current.getAmbigBases());
					pDataArray->longHomoPolymer.push_back(current.getLongHomoPolymer());
				}
            }
		}
		
		in.close();
		
		return 0;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "ScreenSeqsCommand", "MySumThreadFunction");
		exit(1);
	}
} 

/**************************************************************************************************/

static DWORD WINAPI MySumScreenThreadFunction(LPVOID lpParam){ 
	sumScreenData* pDataArray;
	pDataArray = (sumScreenData*)lpParam;
	
	try {
        
        ofstream goodFile;
		pDataArray->m->openOutputFile(pDataArray->goodFName, goodFile);
		
		ofstream badAccnosFile;
		pDataArray->m->openOutputFile(pDataArray->badAccnosFName, badAccnosFile);
		
		ifstream in;
		pDataArray->m->openInputFile(pDataArray->filename, in);
        
		//print header if you are process 0
		if ((pDataArray->start == 0) || (pDataArray->start == 1)) {
			in.seekg(0);
		}else { //this accounts for the difference in line endings. 
			in.seekg(pDataArray->start-1); pDataArray->m->gobble(in); 
		}
		
		pDataArray->count = pDataArray->end;
		for(int i = 0; i < pDataArray->end; i++){ //end is the number of sequences to process
			
			if (pDataArray->m->control_pressed) { in.close(); badAccnosFile.close(); goodFile.close(); pDataArray->count = 1; return 1; }
			
			Sequence currSeq(in); pDataArray->m->gobble(in); 
			
			if (currSeq.getName() != "") {
				bool goodSeq = 1;		//	innocent until proven guilty
				if(goodSeq == 1 && pDataArray->startPos != -1 && pDataArray->startPos < currSeq.getStartPos())			{	goodSeq = 0;	}
				if(goodSeq == 1 && pDataArray->endPos != -1 && pDataArray->endPos > currSeq.getEndPos())				{	goodSeq = 0;	}
				if(goodSeq == 1 && pDataArray->maxAmbig != -1 && pDataArray->maxAmbig <	currSeq.getAmbigBases())		{	goodSeq = 0;	}
				if(goodSeq == 1 && pDataArray->maxHomoP != -1 && pDataArray->maxHomoP < currSeq.getLongHomoPolymer())	{	goodSeq = 0;	}
				if(goodSeq == 1 && pDataArray->minLength != -1 && pDataArray->minLength > currSeq.getNumBases())		{	goodSeq = 0;	}
				if(goodSeq == 1 && pDataArray->maxLength != -1 && pDataArray->maxLength < currSeq.getNumBases())		{	goodSeq = 0;	}
				
				if(goodSeq == 1){
					currSeq.printSequence(goodFile);	
				}
				else{
					badAccnosFile << currSeq.getName() << endl;
					pDataArray->badSeqNames->insert(currSeq.getName());
				}
    
			}		
            //report progress
			if((i+1) % 100 == 0){	pDataArray->m->mothurOut("Processing sequence: " + toString(i+1)); pDataArray->m->mothurOutEndLine();		}
		}
		//report progress
		if((pDataArray->count) % 100 != 0){	pDataArray->m->mothurOut("Processing sequence: " + toString(pDataArray->count)); pDataArray->m->mothurOutEndLine();		}
		

		
		in.close();
        goodFile.close();
        badAccnosFile.close();
		
		return 0;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "ScreenSeqsCommand", "MySumScreenThreadFunction");
		exit(1);
	}
} 

#endif

/**************************************************************************************************/



#endif
