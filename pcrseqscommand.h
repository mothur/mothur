#ifndef Mothur_pcrseqscommand_h
#define Mothur_pcrseqscommand_h

//
//  pcrseqscommand.h
//  Mothur
//
//  Created by Sarah Westcott on 3/14/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//


#include "command.hpp"
#include "sequence.hpp"
#include "trimoligos.h"
#include "alignment.hpp"
#include "needlemanoverlap.hpp"
#include "counttable.h"

class PcrSeqsCommand : public Command {
public:
	PcrSeqsCommand(string);
	PcrSeqsCommand();
	~PcrSeqsCommand(){}
	
	vector<string> setParameters();
	string getCommandName()			{ return "pcr.seqs";	}
	string getCommandCategory()		{ return "Sequence Processing";		}
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "http://www.mothur.org/wiki/Pcr.seqs"; }
	string getDescription()		{ return "pcr.seqs"; }
    
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
    
    struct linePair {
        unsigned long long start;
        unsigned long long end;
        linePair(unsigned long long i, unsigned long long j) : start(i), end(j) {}
        linePair() {}
    };
    
    vector<linePair> lines;
	bool getOligos(vector<vector<string> >&, vector<vector<string> >&, vector<vector<string> >&);
    bool abort, keepprimer, keepdots;
	string fastafile, oligosfile, taxfile, groupfile, namefile, countfile, ecolifile, outputDir, nomatch;
	int start, end, processors, length, pdiffs;
	
    vector<string> revPrimer, outputNames;
	map<string, int> primers;
    
    int writeAccnos(set<string>);
    int readName(set<string>&);
    int readGroup(set<string>);
    int readTax(set<string>);
    int readCount(set<string>);
    bool readOligos();
    bool readEcoli();
	int driverPcr(string, string, string, set<string>&, linePair);	
	int createProcesses(string, string, string, set<string>&);
    bool isAligned(string, map<int, int>&);
    string reverseOligo(string);
};

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct pcrData {
	string filename; 
    string goodFasta, badFasta, oligosfile, ecolifile, nomatch;
	unsigned long long fstart;
	unsigned long long fend;
	int count, start, end, length, pdiffs;
	MothurOut* m;
	map<string, int> primers;
    vector<string> revPrimer;
    set<string> badSeqNames;
    bool keepprimer, keepdots;
	
	
	pcrData(){}
	pcrData(string f, string gf, string bfn, MothurOut* mout, string ol, string ec, map<string, int> pr, vector<string> rpr, string nm, bool kp, bool kd, int st, int en, int l, int pd, unsigned long long fst, unsigned long long fen) {
		filename = f;
        goodFasta = gf;
        badFasta = bfn;
		m = mout;
        oligosfile = ol;
        ecolifile = ec;
        primers = pr;
        revPrimer = rpr;
        nomatch = nm;
        keepprimer = kp;
        keepdots = kd;
		start = st;
		end = en;
        length = l;
		fstart = fst;
        fend = fen;
        pdiffs = pd;
		count = 0;
	}
};
/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyPcrThreadFunction(LPVOID lpParam){ 
	pcrData* pDataArray;
	pDataArray = (pcrData*)lpParam;
	
	try {
        ofstream goodFile;
		pDataArray->m->openOutputFile(pDataArray->goodFasta, goodFile);
        
        ofstream badFile;
		pDataArray->m->openOutputFile(pDataArray->badFasta, badFile);
		
		ifstream inFASTA;
		pDataArray->m->openInputFile(pDataArray->filename, inFASTA);
        
		//print header if you are process 0
		if ((pDataArray->fstart == 0) || (pDataArray->fstart == 1)) {
			inFASTA.seekg(0);
		}else { //this accounts for the difference in line endings. 
			inFASTA.seekg(pDataArray->fstart-1); pDataArray->m->gobble(inFASTA); 
		}
        
        set<int> lengths;
        //pdiffs, bdiffs, primers, barcodes, revPrimers
        map<string, int> faked;
        TrimOligos trim(pDataArray->pdiffs, 0, pDataArray->primers, faked, pDataArray->revPrimer);
		
		for(int i = 0; i < pDataArray->fend; i++){ //end is the number of sequences to process
            pDataArray->count++;
			if (pDataArray->m->control_pressed) {  break; }
			
			Sequence currSeq(inFASTA); pDataArray->m->gobble(inFASTA);
          
            string trashCode = "";
			if (currSeq.getName() != "") {
                
                bool goodSeq = true;
                if (pDataArray->oligosfile != "") {
                    map<int, int> mapAligned;
                    //bool aligned = isAligned(currSeq.getAligned(), mapAligned);
                    ///////////////////////////////////////////////////////////////
                    bool aligned = false;
                    string seq = currSeq.getAligned(); 
                    int countBases = 0;
                    for (int k = 0; k < seq.length(); k++) {
                        if (!isalpha(seq[k])) { aligned = true; }
                        else { mapAligned[countBases] = k; countBases++; } //maps location in unaligned -> location in aligned.
                    }                                                   //ie. the 3rd base may be at spot 10 in the alignment
                                                                        //later when we trim we want to trim from spot 10.
                    ///////////////////////////////////////////////////////////////
                    
                    //process primers
                    if (pDataArray->primers.size() != 0) {
                        int primerStart = 0; int primerEnd = 0;
                        bool good = trim.findForward(currSeq, primerStart, primerEnd);
                        
                        if(!good){	if (pDataArray->nomatch == "reject") { goodSeq = false; } trashCode += "f";	}
                        else{
                            //are you aligned
                            if (aligned) { 
                                if (!pDataArray->keepprimer)    {  
                                    if (pDataArray->keepdots)   { currSeq.filterToPos(mapAligned[primerEnd]);   }
                                    else            { currSeq.setAligned(currSeq.getAligned().substr(mapAligned[primerEnd]));                                              }
                                } 
                                else                {  
                                    if (pDataArray->keepdots)   { currSeq.filterToPos(mapAligned[primerStart]);  }
                                    else            { currSeq.setAligned(currSeq.getAligned().substr(mapAligned[primerStart]));                                              }
                                }
                                ///////////////////////////////////////////////////////////////
                                mapAligned.clear();
                                string seq = currSeq.getAligned(); 
                                int countBases = 0;
                                for (int k = 0; k < seq.length(); k++) {
                                    if (!isalpha(seq[k])) { ; }
                                    else { mapAligned[countBases] = k; countBases++; } 
                                }                                                   
                                ///////////////////////////////////////////////////////////////
                            }else { 
                                if (!pDataArray->keepprimer)    { currSeq.setAligned(currSeq.getUnaligned().substr(primerEnd)); } 
                                else                { currSeq.setAligned(currSeq.getUnaligned().substr(primerStart)); } 
                            }
                        }
                    }
                    
                    //process reverse primers
                    if (pDataArray->revPrimer.size() != 0) {
                        int primerStart = 0; int primerEnd = 0;
                        bool good = trim.findReverse(currSeq, primerStart, primerEnd);
                         
                        if(!good){	if (pDataArray->nomatch == "reject") { goodSeq = false; } trashCode += "r";	}
                        else{ 
                            //are you aligned
                            if (aligned) { 
                                if (!pDataArray->keepprimer)    {  
                                    if (pDataArray->keepdots)   { currSeq.filterFromPos(mapAligned[primerStart]); }
                                    else            { currSeq.setAligned(currSeq.getAligned().substr(0, mapAligned[primerStart]));   }
                                } 
                                else                {  
                                    if (pDataArray->keepdots)   { currSeq.filterFromPos(mapAligned[primerEnd]); }
                                    else            { currSeq.setAligned(currSeq.getAligned().substr(0, mapAligned[primerEnd]));   }
                                }                             }
                            else { 
                                if (!pDataArray->keepprimer)    { currSeq.setAligned(currSeq.getUnaligned().substr(0, primerStart));   } 
                                else                { currSeq.setAligned(currSeq.getUnaligned().substr(0, primerEnd));     }
                            }
                        }
                    }
                }else if (pDataArray->ecolifile != "") {
                    //make sure the seqs are aligned
                    lengths.insert(currSeq.getAligned().length());
                    if (lengths.size() > 1) { pDataArray->m->mothurOut("[ERROR]: seqs are not aligned. When using start and end your sequences must be aligned.\n"); pDataArray->m->control_pressed = true; break; }
                    else if (currSeq.getAligned().length() != pDataArray->length) {
                        pDataArray->m->mothurOut("[ERROR]: seqs are not the same length as ecoli seq. When using ecoli option your sequences must be aligned and the same length as the ecoli sequence.\n"); pDataArray->m->control_pressed = true; break; 
                    }else {
                        if (pDataArray->keepdots)   { 
                            currSeq.filterToPos(pDataArray->start); 
                            currSeq.filterFromPos(pDataArray->end);
                        }else {
                            string seqString = currSeq.getAligned().substr(0, pDataArray->end);
                            seqString = seqString.substr(pDataArray->start);
                            currSeq.setAligned(seqString); 
                        }
                    }
                }else{ //using start and end to trim
                    //make sure the seqs are aligned
                    lengths.insert(currSeq.getAligned().length());
                    if (lengths.size() > 1) { pDataArray->m->mothurOut("[ERROR]: seqs are not aligned. When using start and end your sequences must be aligned.\n"); pDataArray->m->control_pressed = true; break; }
                    else {
                        if (pDataArray->end != -1) {
                            if (pDataArray->end > currSeq.getAligned().length()) {  pDataArray->m->mothurOut("[ERROR]: end is longer than your sequence length, aborting.\n"); pDataArray->m->control_pressed = true; break; }
                            else {
                                if (pDataArray->keepdots)   { currSeq.filterFromPos(pDataArray->end); }
                                else {
                                    string seqString = currSeq.getAligned().substr(0, pDataArray->end);
                                    currSeq.setAligned(seqString); 
                                }
                            }
                        }
                        if (pDataArray->start != -1) { 
                            if (pDataArray->keepdots)   {  currSeq.filterToPos(pDataArray->start);  }
                            else {
                                string seqString = currSeq.getAligned().substr(pDataArray->start);
                                currSeq.setAligned(seqString); 
                            }
                        }
                        
                    }
                }
                
				if(goodSeq == 1)    {   currSeq.printSequence(goodFile);        }
				else {  
                    pDataArray->badSeqNames.insert(currSeq.getName()); 
                    currSeq.setName(currSeq.getName() + '|' + trashCode);
                    currSeq.printSequence(badFile); 
                }
			}
						
			//report progress
			if((i+1) % 100 == 0){	pDataArray->m->mothurOut("Processing sequence: " + toString(i+1)); pDataArray->m->mothurOutEndLine();		}
		}
		//report progress
		if((pDataArray->count) % 100 != 0){	pDataArray->m->mothurOut("Thread Processing sequence: " + toString(pDataArray->count)); pDataArray->m->mothurOutEndLine();		}
		
		goodFile.close();
		inFASTA.close();
        badFile.close();
        
        return 0;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "PcrSeqsCommand", "MyPcrThreadFunction");
		exit(1);
	}
} 

#endif

/**************************************************************************************************/



#endif
