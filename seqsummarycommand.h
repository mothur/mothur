#ifndef SEQSUMMARYCOMMAND_H
#define SEQSUMMARYCOMMAND_H

/*
 *  seqcoordcommand.h
 *  Mothur
 *
 *  Created by Pat Schloss on 5/30/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "sequence.hpp"

/**************************************************************************************************/

class SeqSummaryCommand : public Command {
public:
	SeqSummaryCommand(string);
	SeqSummaryCommand();
	~SeqSummaryCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "summary.seqs";			}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Summary.seqs"; }
	string getDescription()		{ return "summarize the quality of sequences in an unaligned or aligned fasta file"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }		
private:
	bool abort;
	string fastafile, outputDir, namefile, countfile;
	int processors;
	vector<string> outputNames;
	map<string, int> nameMap;
	
	struct linePair {
		unsigned long long start;
		unsigned long long end;
		linePair(unsigned long long i, unsigned long long j) : start(i), end(j) {}
	};

	vector<linePair*> lines;
	vector<int> processIDS;
	
	long long createProcessesCreateSummary(map<int, long long>&, map<int,  long long>&, map<int,  long long>&, map<int,  long long>&, map<int,  long long>&, string, string);
	long long driverCreateSummary(map<int, long long>&, map<int,  long long>&, map<int,  long long>&, map<int,  long long>&, map<int,  long long>&, string, string, linePair*);


};

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct seqSumData {
	map<int, long long> startPosition;
    map<int, long long> endPosition;
    map<int, long long> seqLength;
    map<int, long long> ambigBases;
    map<int, long long> longHomoPolymer;
	string filename; 
	string sumFile; 
	unsigned long long start;
	unsigned long long end;
	int count;
	MothurOut* m;
	bool hasNameMap;
	map<string, int> nameMap;
	
	
	seqSumData(){}
	seqSumData(string f, string sf, MothurOut* mout, unsigned long long st, unsigned long long en, bool na, map<string, int> nam) {
		filename = f;
		sumFile = sf;
		m = mout;
		start = st;
		end = en;
		hasNameMap = na;
		nameMap = nam;
		count = 0;
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MySeqSumThreadFunction(LPVOID lpParam){ 
	seqSumData* pDataArray;
	pDataArray = (seqSumData*)lpParam;
	
	try {
		ofstream outSummary;
		pDataArray->m->openOutputFile(pDataArray->sumFile, outSummary);
		
		ifstream in;
		pDataArray->m->openInputFile(pDataArray->filename, in);

		//print header if you are process 0
		if ((pDataArray->start == 0) || (pDataArray->start == 1)) {
			outSummary << "seqname\tstart\tend\tnbases\tambigs\tpolymer\tnumSeqs" << endl;	
			in.seekg(0);
		}else { //this accounts for the difference in line endings. 
			in.seekg(pDataArray->start-1); pDataArray->m->gobble(in); 
		}
		
		for(int i = 0; i < pDataArray->end; i++){ //end is the number of sequences to process
			
            pDataArray->count++;
            
			if (pDataArray->m->control_pressed) { in.close(); outSummary.close(); pDataArray->count = 1; return 1; }
			
			Sequence current(in); pDataArray->m->gobble(in); 
			
			if (current.getName() != "") {
				
				int num = 1;
				if (pDataArray->hasNameMap){
					//make sure this sequence is in the namefile, else error 
					map<string, int>::iterator it = pDataArray->nameMap.find(current.getName());
					
					if (it == pDataArray->nameMap.end()) { pDataArray->m->mothurOut("[ERROR]: " + current.getName() + " is not in your name or count file, please correct."); pDataArray->m->mothurOutEndLine(); pDataArray->m->control_pressed = true; }
					else { num = it->second; }
				}
				
				int thisStartPosition = current.getStartPos();
                map<int, long long>::iterator it = pDataArray->startPosition.find(thisStartPosition);
                if (it == pDataArray->startPosition.end()) { pDataArray->startPosition[thisStartPosition] = num; } //first finding of this start position, set count.
                else { it->second += num; } //add counts
                
                int thisEndPosition = current.getEndPos();
                it = pDataArray->endPosition.find(thisEndPosition);
                if (it == pDataArray->endPosition.end()) { pDataArray->endPosition[thisEndPosition] = num; } //first finding of this end position, set count.
                else { it->second += num; } //add counts
                
                int thisSeqLength = current.getNumBases();
                it = pDataArray->seqLength.find(thisSeqLength);
                if (it == pDataArray->seqLength.end()) { pDataArray->seqLength[thisSeqLength] = num; } //first finding of this length, set count.
                else { it->second += num; } //add counts
                
                int thisAmbig = current.getAmbigBases();
                it = pDataArray->ambigBases.find(thisAmbig);
                if (it == pDataArray->ambigBases.end()) { pDataArray->ambigBases[thisAmbig] = num; } //first finding of this ambig, set count.
                else { it->second += num; } //add counts
                
                int thisHomoP = current.getLongHomoPolymer();
                it = pDataArray->longHomoPolymer.find(thisHomoP);
                if (it == pDataArray->longHomoPolymer.end()) { pDataArray->longHomoPolymer[thisHomoP] = num; } //first finding of this homop, set count.
                else { it->second += num; } //add counts
                
				pDataArray->count++;
				outSummary << current.getName() << '\t';
				outSummary << thisStartPosition << '\t' << thisEndPosition << '\t';
				outSummary << thisSeqLength << '\t' << thisAmbig << '\t';
				outSummary << thisHomoP << '\t' << num << endl;
			}
		}
		
		in.close();
		outSummary.close();
		
		return 0;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "SeqSummaryCommand", "MySeqSumThreadFunction");
		exit(1);
	}
} 
#endif




#endif

/**************************************************************************************************/


