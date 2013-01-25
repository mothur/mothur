//
//  primerdesigncommand.h
//  Mothur
//
//  Created by Sarah Westcott on 1/18/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#ifndef Mothur_primerdesigncommand_h
#define Mothur_primerdesigncommand_h

#include "command.hpp"
#include "listvector.hpp"
#include "inputdata.h"
#include "sequence.hpp"
#include "alignment.hpp"
#include "needlemanoverlap.hpp"

/**************************************************************************************************/

class PrimerDesignCommand : public Command {
public:
    PrimerDesignCommand(string);
    PrimerDesignCommand();
    ~PrimerDesignCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "primer.design";		}
    string getCommandCategory()		{ return "OTU-Based Approaches";		} 
    
    string getOutputPattern(string);
	string getHelpString();	
    string getCitation() { return "http://www.mothur.org/wiki/Primer.design"; }
    string getDescription()		{ return "design sequence fragments that are specific to particular OTUs"; }
    
    int execute(); 
    void help() { m->mothurOut(getHelpString()); }	
    
private:
    
    struct linePair {
		int start;
		int end;
		linePair(int i, int j) : start(i), end(j) {}
	};
    struct fastaLinePair {
		unsigned long long start;
		unsigned long long end;
		fastaLinePair(unsigned long long i, unsigned long long j) : start(i), end(j) {}
	};
    
    bool abort, allLines, large;
    int cutoff, pdiffs, length, otunumber, processors, alignedLength;
    string outputDir, listfile, namefile, countfile, fastafile, label;
    double minTM, maxTM;
    ListVector* list;
    vector<string> outputNames;

    int initializeCounts(vector< vector< vector<unsigned int> > >& counts, int length, map<string, int>&, map<string, int>&, vector<unsigned int>&);
    map<string, int> readCount(unsigned long int&);
    char getBase(vector<unsigned int> counts, int size);
    int getListVector();
    int countDiffs(string, string);
    set<string> getPrimer(Sequence);
    bool findPrimer(string, string, vector<int>&, vector<int>&, vector<int>&);
    int findMeltingPoint(string primer, double&, double&);
    
    set<int> createProcesses(string, vector<double>&, vector<double>&, set<string>&, vector<Sequence>&);
    set<int> driver(string, vector<double>&, vector<double>&, set<string>&, vector<Sequence>&, int, int, int&);
    vector< vector< vector<unsigned int> > > driverGetCounts(map<string, int>&, unsigned long int&, vector<unsigned int>&, unsigned long long&, unsigned long long&);
    vector<Sequence> createProcessesConSeqs(map<string, int>&, unsigned long int&);
    
};

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct primerDesignData {
	string summaryFileName;
	MothurOut* m;
	int start;
	int end;
	int pdiffs, threadID, otunumber, length;
	set<string> primers;
	vector<double> minTms, maxTms;
    set<int> otusToRemove;
    vector<Sequence> consSeqs;
    int numBinsProcessed;
	
	primerDesignData(){}
	primerDesignData(string sf, MothurOut* mout, int st, int en, vector<double> min, vector<double> max, set<string> pri, vector<Sequence> seqs, int d, int otun, int l, int tid) {
		summaryFileName = sf;
		m = mout;
		start = st;
		end = en;
		pdiffs = d;
        minTms = min;
        maxTms = max;
        primers = pri;
        consSeqs = seqs;
        otunumber = otun;
        length = l;
		threadID = tid;
        numBinsProcessed = 0;
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyPrimerThreadFunction(LPVOID lpParam){ 
	primerDesignData* pDataArray;
	pDataArray = (primerDesignData*)lpParam;
	
	try {
		ofstream outSum;
        pDataArray->m->openOutputFileAppend(pDataArray->summaryFileName, outSum);
        
        for (int i = pDataArray->start; i < pDataArray->end; i++) {
            
            if (pDataArray->m->control_pressed) { break; }
            
            if (i != (pDataArray->otunumber-1)) {
                int primerIndex = 0;
                for (set<string>::iterator it = pDataArray->primers.begin(); it != pDataArray->primers.end(); it++) {
                    vector<int> primerStarts;
                    vector<int> primerEnds;
                    vector<int> mismatches;
                    
                    //bool found = findPrimer(conSeqs[i].getUnaligned(), (*it), primerStarts, primerEnds, mismatches);
                    ///////////////////////////////////////////////////////////////////////////////////////////////////
                    bool found = false;  //innocent til proven guilty
                    
                    string rawSequence = pDataArray->consSeqs[i].getUnaligned();
                    string primer = *it;
                    
                    //look for exact match
                    if(rawSequence.length() < primer.length()) {  found = false;  }
                    else {
                        //search for primer
                        for (int j = 0; j < rawSequence.length()-pDataArray->length; j++){
                            
                            if (pDataArray->m->control_pressed) {  found = false; break; }
                            
                            string rawChunk = rawSequence.substr(j, pDataArray->length);
                            
                            //int numDiff = countDiffs(primer, rawchuck);
                            ///////////////////////////////////////////////////////////////////////
                            int numDiff = 0;
                            string oligo = primer;
                            string seq = rawChunk;
                            
                            for(int k=0;k<pDataArray->length;k++){
                                
                                oligo[k] = toupper(oligo[k]);
                                seq[k] = toupper(seq[k]);
                               
                                if(oligo[k] != seq[k]){
            
                                    if((oligo[k] == 'N' || oligo[k] == 'I') && (seq[k] == 'N'))				{	numDiff++;	}
                                    else if(oligo[k] == 'R' && (seq[k] != 'A' && seq[k] != 'G'))					{	numDiff++;	}
                                    else if(oligo[k] == 'Y' && (seq[k] != 'C' && seq[k] != 'T'))					{	numDiff++;	}
                                    else if(oligo[k] == 'M' && (seq[k] != 'C' && seq[k] != 'A'))					{	numDiff++;	}
                                    else if(oligo[k] == 'K' && (seq[k] != 'T' && seq[k] != 'G'))					{	numDiff++;	}
                                    else if(oligo[k] == 'W' && (seq[k] != 'T' && seq[k] != 'A'))					{	numDiff++;	}
                                    else if(oligo[k] == 'S' && (seq[k] != 'C' && seq[k] != 'G'))					{	numDiff++;	}
                                    else if(oligo[k] == 'B' && (seq[k] != 'C' && seq[k] != 'T' && seq[k] != 'G'))	{	numDiff++;	}
                                    else if(oligo[k] == 'D' && (seq[k] != 'A' && seq[k] != 'T' && seq[k] != 'G'))	{	numDiff++;	}
                                    else if(oligo[k] == 'H' && (seq[k] != 'A' && seq[k] != 'T' && seq[k] != 'C'))	{	numDiff++;	}
                                    else if(oligo[k] == 'V' && (seq[k] != 'A' && seq[k] != 'C' && seq[k] != 'G'))	{	numDiff++;	}
                                    else if(oligo[k] == 'A' && (seq[k] != 'A' && seq[k] != 'M' && seq[k] != 'R' && seq[k] != 'W' && seq[k] != 'D' && seq[k] != 'H' && seq[k] != 'V'))       {	numDiff++;	}
                                    else if(oligo[k] == 'C' && (seq[k] != 'C' && seq[k] != 'Y' && seq[k] != 'M' && seq[k] != 'S' && seq[k] != 'B' && seq[k] != 'H' && seq[k] != 'V'))       {	numDiff++;	}
                                    else if(oligo[k] == 'G' && (seq[k] != 'G' && seq[k] != 'R' && seq[k] != 'K' && seq[k] != 'S' && seq[k] != 'B' && seq[k] != 'D' && seq[k] != 'V'))       {	numDiff++;	}
                                    else if(oligo[k] == 'T' && (seq[k] != 'T' && seq[k] != 'Y' && seq[k] != 'K' && seq[k] != 'W' && seq[k] != 'B' && seq[k] != 'D' && seq[k] != 'H'))       {	numDiff++;	}
                                    else if((oligo[k] == '.' || oligo[k] == '-'))           {	numDiff++;	}
                                }
                            }
                            ///////////////////////////////////////////////////////////////////////
                            
                            if(numDiff <= pDataArray->pdiffs){
                                primerStarts.push_back(j);
                                primerEnds.push_back(j+pDataArray->length);
                                mismatches.push_back(numDiff);
                                found = true;
                            }
                        }
                    }
                    ///////////////////////////////////////////////////////////////////////////////////////////////////
                    
                    //if we found it report to the table
                    if (found) {
                        for (int j = 0; j < primerStarts.size(); j++) {
                            outSum << (i+1) << '\t' << *it << '\t' << primerStarts[j] << '\t' << primerEnds[j] << '\t' << pDataArray->length << '\t' << mismatches[j] << '\t' << pDataArray->minTms[primerIndex] << '\t' << pDataArray->maxTms[primerIndex] << endl;
                        }
                        pDataArray->otusToRemove.insert(i);
                    }
                    primerIndex++;
                }
            }
            pDataArray->numBinsProcessed++;
        }
        outSum.close();
        
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "PrimerDesignCommand", "MyPrimerThreadFunction");
		exit(1);
	}
} 
#endif

/**************************************************************************************************/





#endif
