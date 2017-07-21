#ifndef SUMMARYQUALCOMMAND_H
#define SUMMARYQUALCOMMAND_H

/*
 *  summaryqualcommand.h
 *  Mothur
 *
 *  Created by westcott on 11/28/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"
#include "qualityscores.h"

/**************************************************************************************************/

class SummaryQualCommand : public Command {
public:
	SummaryQualCommand(string);
	SummaryQualCommand();
	~SummaryQualCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "summary.qual";			}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Summary.qual"; }
	string getDescription()		{ return "summarize the quality of a set of sequences"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	bool abort;
	string qualfile, outputDir, namefile, countfile;
	vector<string> outputNames;
	map<string, int> nameMap;
	int processors;
	
	vector<linePair> lines;
	vector<int> processIDS;
	
	int createProcessesCreateSummary(vector<int>&, vector<int>&, vector< vector<int> >&, string);
	int driverCreateSummary(vector<int>&, vector<int>&, vector< vector<int> >&, string, linePair);	
	int printQual(string, vector<int>&, vector<int>&, vector< vector<int> >&);
};

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct seqSumQualData {
	vector<int> position;
	vector<int> averageQ;
	vector< vector<int> > scores; 
	string filename; 
	unsigned long long start;
	unsigned long long end;
	int count, numSeqs;
	MothurOut* m;
    bool hasNameMap;
	map<string, int> nameMap;
	
	~seqSumQualData(){}
	seqSumQualData(string f, MothurOut* mout, unsigned long long st, unsigned long long en, bool n, map<string, int> nam) {
		filename = f;
		m = mout;
		start = st;
		end = en;
		hasNameMap = n;
		nameMap = nam;
		count = 0;
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MySeqSumQualThreadFunction(LPVOID lpParam){ 
	seqSumQualData* pDataArray;
	pDataArray = (seqSumQualData*)lpParam;
	
	try {
		ifstream in;
		pDataArray->m->openInputFile(pDataArray->filename, in);
		
		//print header if you are process 0
		if ((pDataArray->start == 0) || (pDataArray->start == 1)) {
			in.seekg(0);
            pDataArray->m->zapGremlins(in);
		}else { //this accounts for the difference in line endings. 
			in.seekg(pDataArray->start-1); pDataArray->m->gobble(in); 
		}
		
		pDataArray->count = 0;
        pDataArray->numSeqs = 0;
		for(int i = 0; i < pDataArray->end; i++){ //end is the number of sequences to process
				
			if (pDataArray->m->control_pressed) { in.close(); pDataArray->count = 1; return 1; }
			
			QualityScores current(in); pDataArray->m->gobble(in);
			
			if (current.getName() != "") {
			
				int num = 1;
				if (pDataArray->hasNameMap) {
					//make sure this sequence is in the namefile, else error 
					map<string, int>::iterator it = pDataArray->nameMap.find(current.getName());
					
					if (it == pDataArray->nameMap.end()) { pDataArray->m->mothurOut("[ERROR]: " + current.getName() + " is not in your namefile, please correct."); pDataArray->m->mothurOutEndLine(); pDataArray->m->control_pressed = true; }
					else { num = it->second; }
				}
				
				vector<int> thisScores = current.getQualityScores();
				
				//resize to num of positions setting number of seqs with that size to 1
				if (pDataArray->position.size() < thisScores.size()) { pDataArray->position.resize(thisScores.size(), 0); }
				if (pDataArray->averageQ.size() < thisScores.size()) { pDataArray->averageQ.resize(thisScores.size(), 0); }
				if (pDataArray->scores.size() < thisScores.size()) { 
					pDataArray->scores.resize(thisScores.size()); 
					for (int i = 0; i < pDataArray->scores.size(); i++) { pDataArray->scores.at(i).resize(41, 0); }
				}
				
				//increase counts of number of seqs with this position
				//average is really the total, we will average in execute
				for (int i = 0; i < thisScores.size(); i++) { 
					pDataArray->position.at(i) += num; 
					pDataArray->averageQ.at(i) += (thisScores[i] * num); //weighting for namesfile
					if (thisScores[i] > 41) { pDataArray->m->mothurOut("[ERROR]: " + current.getName() + " has a quality scores of " + toString(thisScores[i]) + ", expecting values to be less than 40."); pDataArray->m->mothurOutEndLine(); pDataArray->m->control_pressed = true; }
					else { pDataArray->scores.at(i)[thisScores[i]] += num; }  
				}
				
				pDataArray->numSeqs += num;
                pDataArray->count++;
			}
		}
		
		in.close();
		
		return 0;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "SummaryQualCommand", "MySeqSumQualThreadFunction");
		exit(1);
	}
} 
#endif


/**************************************************************************************************/


#endif

