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
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Screen.seqs"; }
	string getDescription()		{ return "enables you to keep sequences that fulfill certain user defined criteria"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
	
private:
	vector<linePair> lines;

    int optimizeContigs();
    int optimizeAlign();
	int driver(linePair, string, string, string, map<string, string>&);
	int createProcesses(string, string, string, map<string, string>&);
    int screenSummary(map<string, string>&);
    int screenContigs(map<string, string>&);
    int screenAlignReport(map<string, string>&);
    int runFastaScreening(map<string, string>&);
    int screenFasta(map<string, string>&);
    int screenReports(map<string, string>&);
	int getSummary();
    int getSummaryReport();
 
    bool abort;
    string fastafile, namefile, groupfile, alignreport, outputDir, qualfile, taxonomy, countfile, contigsreport, summaryfile, fileType, badAccnosFile;
	int startPos, endPos, maxAmbig, maxHomoP, minLength, maxLength, processors, minOverlap, oStart, oEnd, mismatches, maxN, maxInsert;
    float minSim, minScore, criteria;
	vector<string> outputNames;
	vector<string> optimize;
	map<string, int> nameMap;
	
    
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
    vector<int> numNs;
	string filename, namefile, countfile; 
	unsigned long long start;
	unsigned long long end;
	int count;
	MothurOut* m;
	map<string, int> nameMap;
	
	
	sumData(){}
	sumData(string f, MothurOut* mout, unsigned long long st, unsigned long long en, string nf, string cf, map<string, int> nam) {
		filename = f;
        namefile = nf;
        countfile = cf;
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
    int startPos, endPos, maxAmbig, maxHomoP, minLength, maxLength, maxN;
	unsigned long long start;
	unsigned long long end;
	int count;
	MothurOut* m;
	string goodFName, badAccnosFName, filename;
    map<string, string> badSeqNames;
    string summaryfile, contigsreport;
	
	
	sumScreenData(){}
	sumScreenData(int s, int e, int a, int h, int minl, int maxl, int mn, map<string, string> bs, string f, string sum, string cont, MothurOut* mout, unsigned long long st, unsigned long long en, string gf, string bf) {
		startPos = s;
		endPos = e;
		minLength = minl;
        maxLength = maxl;
		maxAmbig = a;
		maxHomoP = h;
        maxN = mn;
		filename = f;
        goodFName = gf;
        badAccnosFName = bf;
		m = mout;
		start = st;
		end = en;
        summaryfile = sum;
        contigsreport = cont;
        badSeqNames = bs;
		count = 0;
	}
};
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
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
            pDataArray->m->zapGremlins(in);
		}else { //this accounts for the difference in line endings. 
			in.seekg(pDataArray->start-1); pDataArray->m->gobble(in); 
		}
        
		for(int i = 0; i < pDataArray->end; i++){ //end is the number of sequences to process
			
            pDataArray->count++;
            
			if (pDataArray->m->control_pressed) { in.close(); badAccnosFile.close(); goodFile.close(); pDataArray->count = 1; return 1; }
			
			Sequence currSeq(in); pDataArray->m->gobble(in); 
			
			if (currSeq.getName() != "") {
				bool goodSeq = 1;		//	innocent until proven guilty
                string trashCode = "";
                //have the report files found you bad
                map<string, string>::iterator it = pDataArray->badSeqNames.find(currSeq.getName());
                if (it != pDataArray->badSeqNames.end()) { goodSeq = 0;  trashCode = it->second; } //found it 
                
                if (pDataArray->summaryfile == "") {
                    if(pDataArray->startPos != -1 && pDataArray->startPos < currSeq.getStartPos())			{	goodSeq = 0;	trashCode += "start|"; }
                    if(pDataArray->endPos != -1 && pDataArray->endPos > currSeq.getEndPos())				{	goodSeq = 0;	trashCode += "end|"; }
                    if(pDataArray->maxAmbig != -1 && pDataArray->maxAmbig <	currSeq.getAmbigBases())		{	goodSeq = 0;	trashCode += "ambig|"; }
                    if(pDataArray->maxHomoP != -1 && pDataArray->maxHomoP < currSeq.getLongHomoPolymer())	{	goodSeq = 0;	trashCode += "homop|"; }
                    if(pDataArray->minLength > currSeq.getNumBases())		{	goodSeq = 0;	trashCode += "<length|"; }
                    if(pDataArray->maxLength != -1 && pDataArray->maxLength < currSeq.getNumBases())		{	goodSeq = 0;	trashCode += ">length|"; }
                }
                if (pDataArray->contigsreport == "") { //contigs report includes this so no need to check again
                    if(pDataArray->maxN != -1 && pDataArray->maxN < currSeq.getNumNs())                     {	goodSeq = 0;	trashCode += "n|"; }
                }
				
                
				if(goodSeq == 1){
					currSeq.printSequence(goodFile);	
				}
				else{
					badAccnosFile << currSeq.getName() << '\t' << trashCode.substr(0, trashCode.length()-1) << endl;
					pDataArray->badSeqNames[currSeq.getName()] = trashCode;
				}
    
			}		
            //report progress
			if((i+1) % 100 == 0){	pDataArray->m->mothurOutJustToScreen("Processing sequence: " + toString(i+1)+"\n"); 		}
		}
		//report progress
		if((pDataArray->count) % 100 != 0){	pDataArray->m->mothurOutJustToScreen("Processing sequence: " + toString(pDataArray->count)+"\n"); 		}
		

		
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
