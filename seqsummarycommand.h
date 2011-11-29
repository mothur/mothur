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
	string getCitation() { return "http://www.mothur.org/wiki/Summary.seqs"; }
	string getDescription()		{ return "summarize the quality of sequences in an unaligned or aligned fasta file"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }		
private:
	bool abort;
	string fastafile, outputDir, namefile;
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
	
	int createProcessesCreateSummary(vector<int>&, vector<int>&, vector<int>&, vector<int>&, vector<int>&, string, string);
	int driverCreateSummary(vector<int>&, vector<int>&, vector<int>&, vector<int>&, vector<int>&, string, string, linePair*);	

	#ifdef USE_MPI
	int MPICreateSummary(int, int, vector<int>&, vector<int>&, vector<int>&, vector<int>&, vector<int>&, MPI_File&, MPI_File&, vector<unsigned long long>&);	
	#endif


};

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct seqSumData {
	vector<int>* startPosition;
	vector<int>* endPosition;
	vector<int>* seqLength; 
	vector<int>* ambigBases; 
	vector<int>* longHomoPolymer; 
	string filename; 
	string sumFile; 
	unsigned long long start;
	unsigned long long end;
	int count;
	MothurOut* m;
	string namefile;
	map<string, int> nameMap;
	
	
	seqSumData(){}
	seqSumData(vector<int>* s, vector<int>* e, vector<int>* l, vector<int>* a, vector<int>* h, string f, string sf, MothurOut* mout, unsigned long long st, unsigned long long en, string na, map<string, int> nam) {
		startPosition = s;
		endPosition = e;
		seqLength = l;
		ambigBases = a;
		longHomoPolymer = h;
		filename = f;
		sumFile = sf;
		m = mout;
		start = st;
		end = en;
		namefile = na;
		nameMap = nam;
		count = 0;
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
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
		
		pDataArray->count = pDataArray->end;
		for(int i = 0; i < pDataArray->end; i++){ //end is the number of sequences to process
			
			if (pDataArray->m->control_pressed) { in.close(); outSummary.close(); pDataArray->count = 1; return 1; }
			
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
					pDataArray->startPosition->push_back(current.getStartPos());
					pDataArray->endPosition->push_back(current.getEndPos());
					pDataArray->seqLength->push_back(current.getNumBases());
					pDataArray->ambigBases->push_back(current.getAmbigBases());
					pDataArray->longHomoPolymer->push_back(current.getLongHomoPolymer());
				}
				
				outSummary << current.getName() << '\t';
				outSummary << current.getStartPos() << '\t' << current.getEndPos() << '\t';
				outSummary << current.getNumBases() << '\t' << current.getAmbigBases() << '\t';
				outSummary << current.getLongHomoPolymer() << '\t' << num << endl;
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


