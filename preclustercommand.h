#ifndef PRECLUSTERCOMMAND_H
#define PRECLUSTERCOMMAND_H


/*
 *  preclustercommand.h
 *  Mothur
 *
 *  Created by westcott on 12/21/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"
#include "sequence.hpp"
#include "sequenceparser.h"

/************************************************************/
struct seqPNode {
	int numIdentical;
	Sequence seq;
	string names;
	bool active;
	seqPNode() {}
	seqPNode(int n, Sequence s, string nm) : numIdentical(n), seq(s), names(nm), active(1) {}
	~seqPNode() {}
};
/************************************************************/
inline bool comparePriority(seqPNode first, seqPNode second) {  return (first.numIdentical > second.numIdentical); }
//************************************************************/

class PreClusterCommand : public Command {
	
public:
	PreClusterCommand(string);
	PreClusterCommand();
	~PreClusterCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "pre.cluster";				}
	string getCommandCategory()		{ return "Sequence Processing";		}
	string getHelpString();	
	string getCitation() { return "http://www.mothur.org/wiki/Pre.cluster"; }
	string getDescription()		{ return "implements a pseudo-single linkage algorithm with the goal of removing sequences that are likely due to pyrosequencing errors"; }

	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	
	struct linePair {
		int start;
		int end;
		linePair(int i, int j) : start(i), end(j) {}
	};
	
	int diffs, length, processors;
	bool abort, bygroup;
	string fastafile, namefile, outputDir, groupfile;
	vector<seqPNode> alignSeqs; //maps the number of identical seqs to a sequence
	map<string, string> names; //represents the names file first column maps to second column
	map<string, int> sizes;  //this map a seq name to the number of identical seqs in the names file
	map<string, int>::iterator itSize; 
//	map<string, bool> active; //maps sequence name to whether it has already been merged or not.
	vector<string> outputNames;
	map<string, vector<string> > outputTypes;
	
	int readFASTA();
	void readNameFile();
	//int readNamesFASTA();
	int calcMisMatches(string, string);
	void printData(string, string); //fasta filename, names file name
	int process();
	int loadSeqs(map<string, string>&, vector<Sequence>&);
	int driverGroups(SequenceParser*, string, string, int, int, vector<string> groups);
	int createProcessesGroups(SequenceParser*, string, string, vector<string>);
};

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
typedef struct preClusterData {
	string fastafile; 
	string namefile; 
	string groupfile;
	string newFName, newNName;
	MothurOut* m;
	int start;
	int end;
	int diffs, threadID;
	vector<string> groups;
	
	preClusterData(){}
	preClusterData(string f, string n, string g, string nff,  string nnf, vector<string> gr, MothurOut* mout, int st, int en, int d, int tid) {
		fastafile = f;
		namefile = n;
		groupfile = g;
		newFName = nff;
		newNName = nnf;
		m = mout;
		start = st;
		end = en;
		diffs = d;
		threadID = tid;
		groups = gr;
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
#else
static DWORD WINAPI MyPreclusterThreadFunction(LPVOID lpParam){ 
	preClusterData* pDataArray;
	pDataArray = (preClusterData*)lpParam;
	
	try {
		
		//parse fasta and name file by group
		SequenceParser* parser;
		if (pDataArray->namefile != "") { parser = new SequenceParser(pDataArray->groupfile, pDataArray->fastafile, pDataArray->namefile);	}
		else							{ parser = new SequenceParser(pDataArray->groupfile, pDataArray->fastafile);						}
		
		int numSeqs = 0;
		vector<seqPNode> alignSeqs;
		//clear out old files
		ofstream outF; pDataArray->m->openOutputFile(pDataArray->newFName, outF); outF.close();
		ofstream outN; pDataArray->m->openOutputFile(pDataArray->newNName, outN);  outN.close();
		
		//precluster each group
		for (int k = pDataArray->start; k < pDataArray->end; k++) {
			
			int start = time(NULL);
			
			if (pDataArray->m->control_pressed) {  delete parser; return 0; }
			
			pDataArray->m->mothurOutEndLine(); pDataArray->m->mothurOut("Processing group " + pDataArray->groups[k] + ":"); pDataArray->m->mothurOutEndLine();
			
			map<string, string> thisNameMap;
			if (pDataArray->namefile != "") { thisNameMap = parser->getNameMap(pDataArray->groups[k]); }
			vector<Sequence> thisSeqs = parser->getSeqs(pDataArray->groups[k]);
			
			//fill alignSeqs with this groups info.
			////////////////////////////////////////////////////
			//numSeqs = loadSeqs(thisNameMap, thisSeqs); same function below
			
			int length = 0;
			alignSeqs.clear();
			map<string, string>::iterator it;
			bool error = false;
		 	
			for (int i = 0; i < thisSeqs.size(); i++) {
				
				if (pDataArray->m->control_pressed) { delete parser; return 0; }
				
				if (pDataArray->namefile != "") {
					it = thisNameMap.find(thisSeqs[i].getName());
					
					//should never be true since parser checks for this
					if (it == thisNameMap.end()) { pDataArray->m->mothurOut(thisSeqs[i].getName() + " is not in your names file, please correct."); pDataArray->m->mothurOutEndLine(); error = true; }
					else{
						//get number of reps
						int numReps = 1;
						for(int j=0;j<(it->second).length();j++){
							if((it->second)[j] == ','){	numReps++;	}
						}
						
						seqPNode tempNode(numReps, thisSeqs[i], it->second);
						alignSeqs.push_back(tempNode);
						if (thisSeqs[i].getAligned().length() > length) {  length = thisSeqs[i].getAligned().length();  }
					}	
				}else { //no names file, you are identical to yourself 
					seqPNode tempNode(1, thisSeqs[i], thisSeqs[i].getName());
					alignSeqs.push_back(tempNode);
					if (thisSeqs[i].getAligned().length() > length) {  length = thisSeqs[i].getAligned().length();  }
				}
			}
			
			//sanity check
			if (error) { pDataArray->m->control_pressed = true; }
			
			thisSeqs.clear();
			numSeqs = alignSeqs.size();
			
			////////////////////////////////////////////////////
			
			if (pDataArray->m->control_pressed) {   delete parser; return 0; }
			
			if (pDataArray->diffs > length) { pDataArray->m->mothurOut("Error: diffs is greater than your sequence length."); pDataArray->m->mothurOutEndLine(); pDataArray->m->control_pressed = true; return 0;  }
			
			////////////////////////////////////////////////////
			//int count = process(); - same function below
			
			//sort seqs by number of identical seqs
			sort(alignSeqs.begin(), alignSeqs.end(), comparePriority);
			
			int count = 0;
			
			//think about running through twice...
			for (int i = 0; i < numSeqs; i++) {
				
				//are you active
				//			itActive = active.find(alignSeqs[i].seq.getName());
				
				if (alignSeqs[i].active) {  //this sequence has not been merged yet
					
					//try to merge it with all smaller seqs
					for (int j = i+1; j < numSeqs; j++) {
						
						if (pDataArray->m->control_pressed) { delete parser; return 0; }
						
						if (alignSeqs[j].active) {  //this sequence has not been merged yet
							//are you within "diff" bases
							//int mismatch = calcMisMatches(alignSeqs[i].seq.getAligned(), alignSeqs[j].seq.getAligned());
							int mismatch = 0;
							
							for (int k = 0; k < alignSeqs[i].seq.getAligned().length(); k++) {
								//do they match
								if (alignSeqs[i].seq.getAligned()[k] != alignSeqs[j].seq.getAligned()[k]) { mismatch++; }
								if (mismatch > pDataArray->diffs) { mismatch = length; break; } //to far to cluster
							}
							
							if (mismatch <= pDataArray->diffs) {
								//merge
								alignSeqs[i].names += ',' + alignSeqs[j].names;
								alignSeqs[i].numIdentical += alignSeqs[j].numIdentical;
								
								alignSeqs[j].active = 0;
								alignSeqs[j].numIdentical = 0;
								count++;
							}
						}//end if j active
					}//end if i != j
					
					//remove from active list 
					alignSeqs[i].active = 0;
					
				}//end if active i
				if(i % 100 == 0)	{ pDataArray->m->mothurOut(toString(i) + "\t" + toString(numSeqs - count) + "\t" + toString(count)); pDataArray->m->mothurOutEndLine();	}
			}
			
			if(numSeqs % 100 != 0)	{ pDataArray->m->mothurOut(toString(numSeqs) + "\t" + toString(numSeqs - count) + "\t" + toString(count)); pDataArray->m->mothurOutEndLine();	}	
			////////////////////////////////////////////////////
			
			if (pDataArray->m->control_pressed) {  delete parser; return 0; }
			
			pDataArray->m->mothurOut("Total number of sequences before pre.cluster was " + toString(alignSeqs.size()) + ".");pDataArray-> m->mothurOutEndLine();
			pDataArray->m->mothurOut("pre.cluster removed " + toString(count) + " sequences."); pDataArray->m->mothurOutEndLine(); pDataArray->m->mothurOutEndLine(); 
			
			////////////////////////////////////////////////////
			//printData(pDataArray->newFFile, pDataArray->newNFile); - same as below
			ofstream outFasta;
			ofstream outNames;
			
			pDataArray->m->openOutputFileAppend(pDataArray->newFName, outFasta);
			pDataArray->m->openOutputFileAppend(pDataArray->newNName, outNames);
						
			for (int i = 0; i < alignSeqs.size(); i++) {
				if (alignSeqs[i].numIdentical != 0) {
					alignSeqs[i].seq.printSequence(outFasta); 
					outNames << alignSeqs[i].seq.getName() << '\t' << alignSeqs[i].names << endl;
				}
			}
			
			outFasta.close();
			outNames.close();
			////////////////////////////////////////////////////
			
			pDataArray->m->mothurOut("It took " + toString(time(NULL) - start) + " secs to cluster " + toString(numSeqs) + " sequences."); pDataArray->m->mothurOutEndLine(); 
			
		}
		
		return numSeqs;
		

	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "AlignCommand", "MyPreclusterThreadFunction");
		exit(1);
	}
} 
#endif

/**************************************************************************************************/


#endif


