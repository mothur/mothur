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
#include "oligos.h"

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
    vector<linePair> lines;
    bool abort, keepprimer, keepdots, fileAligned, pairedOligos;
	string fastafile, oligosfile, taxfile, groupfile, namefile, countfile, ecolifile, outputDir, nomatch;
	int start, end, processors, length, pdiffs, rdiffs, numFPrimers, numRPrimers;
    Oligos oligos;
	
    vector<string> outputNames;
    
    int writeAccnos(set<string>);
    int readName(set<string>&);
    int readGroup(set<string>);
    int readTax(set<string>);
    int readCount(set<string>);
    int readOligos();
    bool readEcoli();
	int driverPcr(string, string, string, string, set<string>&, linePair, int&, bool&);
	int createProcesses(string, string, string, set<string>&);
    bool isAligned(string, map<int, int>&);
    int adjustDots(string, string, int, int);
    
};

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct pcrData {
	string filename; 
    string goodFasta, badFasta, oligosfile, ecolifile, nomatch, locationsName;
	unsigned long long fstart;
	unsigned long long fend;
	int count, start, end, length, pdiffs, pstart, pend, rdiffs;
	MothurOut* m;
    set<string> badSeqNames;
    bool keepprimer, keepdots, fileAligned, adjustNeeded;
    Utils util;
	
	pcrData(){}
	pcrData(string f, string gf, string bfn, string loc, MothurOut* mout, string ol, string ec, string nm, bool kp, bool kd, int st, int en, int l, int pd, int rd, unsigned long long fst, unsigned long long fen) {
		filename = f;
        goodFasta = gf;
        badFasta = bfn;
		m = mout;
        oligosfile = ol;
        ecolifile = ec;
        nomatch = nm;
        keepprimer = kp;
        keepdots = kd;
        end = en;
		start = st;
        length = l;
		fstart = fst;
        fend = fen;
        pdiffs = pd;
        rdiffs = rd;
        locationsName = loc;
		count = 0;
        fileAligned = true;
        adjustNeeded = false;
        pstart = -1;
        pend = -1;
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
		pDataArray->util.openOutputFile(pDataArray->goodFasta, goodFile);
        
        ofstream badFile;
		pDataArray->util.openOutputFile(pDataArray->badFasta, badFile);
        
        ofstream locationsFile;
		pDataArray->util.openOutputFile(pDataArray->locationsName, locationsFile);
		
		ifstream inFASTA;
		pDataArray->util.openInputFile(pDataArray->filename, inFASTA);
        
		//print header if you are process 0
		if ((pDataArray->fstart == 0) || (pDataArray->fstart == 1)) {
			inFASTA.seekg(0);
		}else { //this accounts for the difference in line endings. 
			inFASTA.seekg(pDataArray->fstart-1); pDataArray->util.gobble(inFASTA); 
		}
        
        set<int> lengths;
        //pdiffs, bdiffs, primers, barcodes, revPrimers
        map<string, int> faked;
        set<int> locations; //locations = beginning locations
        
        Oligos oligos;
        int numFPrimers, numRPrimers;  numFPrimers = 0; numRPrimers = 0;
        map<string, int> primers;
        map<string, int> barcodes; //not used
        vector<string> revPrimer;

        if (pDataArray->oligosfile != "") {
            oligos.read(pDataArray->oligosfile);
            if (oligos.hasPairedPrimers()) {
                map<int, oligosPair> primerPairs = oligos.getPairedPrimers();
                for (map<int, oligosPair>::iterator it = primerPairs.begin(); it != primerPairs.end(); it++) {
                    primers[(it->second).forward] = it->first;
                    revPrimer.push_back((it->second).reverse);
                }
            }else {
                primers = oligos.getPrimers();
                revPrimer = oligos.getReversePrimers();
            }
            numRPrimers = revPrimer.size();
            numFPrimers = primers.size();
        }
        
        TrimOligos trim(pDataArray->pdiffs, pDataArray->rdiffs, 0, primers, barcodes, revPrimer);
		
		for(int i = 0; i < pDataArray->fend; i++){ //end is the number of sequences to process
            pDataArray->count++;
			if (pDataArray->m->getControl_pressed()) {  break; }
			
			Sequence currSeq(inFASTA); pDataArray->util.gobble(inFASTA);
            
            if (pDataArray->fileAligned) { //assume aligned until proven otherwise
                lengths.insert(currSeq.getAligned().length());
                if (lengths.size() > 1) { pDataArray->fileAligned = false; }
            }
            
            string trashCode = "";
            string locationsString = "";
            int thisPStart = -1;
            int thisPEnd = -1;
            int totalDiffs = 0;
            string commentString = "";

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
                    if (numFPrimers != 0) {
                        int primerStart = 0; int primerEnd = 0;
                        vector<int> results = trim.findForward(currSeq, primerStart, primerEnd);
                        bool good = true;
                        if (results[0] > pDataArray->pdiffs) { good = false; }
                        totalDiffs += results[0];
                        commentString += "fpdiffs=" + toString(results[0]) + "(" + trim.getCodeValue(results[1], pDataArray->pdiffs) + ") ";
                        
                        if(!good){	if (pDataArray->nomatch == "reject") { goodSeq = false; } trashCode += "f";	}
                        else{
                            //are you aligned
                            if (aligned) { 
                                if (!pDataArray->keepprimer)    {  
                                    if (pDataArray->keepdots)   { currSeq.filterToPos(mapAligned[primerEnd-1]+1);   }
                                    else            {
                                        currSeq.setAligned(currSeq.getAligned().substr(mapAligned[primerEnd-1]+1));
                                        if (pDataArray->fileAligned) {
                                            thisPStart = mapAligned[primerEnd-1]+1; //locations.insert(mapAligned[primerEnd-1]+1);
                                            locationsString += currSeq.getName() + "\t" + toString(mapAligned[primerEnd-1]+1) + "\n";
                                        }
}
                                } 
                                else                {  
                                    if (pDataArray->keepdots)   { currSeq.filterToPos(mapAligned[primerStart]);  }
                                    else            {
                                        currSeq.setAligned(currSeq.getAligned().substr(mapAligned[primerStart]));
                                        if (pDataArray->fileAligned) {
                                            thisPStart = mapAligned[primerStart]; //locations.insert(mapAligned[primerStart]);
                                            locationsString += currSeq.getName() + "\t" + toString(mapAligned[primerStart]) + "\n";
                                        }
                                    }
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
                    if (numRPrimers != 0) {
                        int primerStart = 0; int primerEnd = 0;
                        vector<int> results = trim.findReverse(currSeq, primerStart, primerEnd);
                        bool good = true;
                        if (results[0] > pDataArray->rdiffs) { good = false; }
                        totalDiffs += results[0];
                        commentString += "rpdiffs=" + toString(results[0]) + "(" + trim.getCodeValue(results[1], pDataArray->rdiffs) + ") ";
                        
                        if(!good){	if (pDataArray->nomatch == "reject") { goodSeq = false; } trashCode += "r";	}
                        else{ 
                            //are you aligned
                            if (aligned) { 
                                if (!pDataArray->keepprimer)    {  
                                    if (pDataArray->keepdots)   { currSeq.filterFromPos(mapAligned[primerStart]); }
                                    else            {
                                        currSeq.setAligned(currSeq.getAligned().substr(0, mapAligned[primerStart]));
                                        if (pDataArray->fileAligned) {
                                            thisPEnd = mapAligned[primerStart]; //locations.insert(mapAligned[primerStart]);
                                            locationsString += currSeq.getName() + "\t" + toString(mapAligned[primerStart]) + "\n";
                                        }

                                    }
                                } 
                                else                {  
                                    if (pDataArray->keepdots)   { currSeq.filterFromPos(mapAligned[primerEnd-1]+1); }
                                    else            {
                                        currSeq.setAligned(currSeq.getAligned().substr(0, mapAligned[primerEnd-1]+1));
                                        if (pDataArray->fileAligned) {
                                            thisPEnd = mapAligned[primerEnd-1]+1; //locations.insert(mapAligned[primerEnd-1]+1);
                                            locationsString += currSeq.getName() + "\t" + toString(mapAligned[primerEnd-1]+1) + "\n";
                                        }

                                    }
                                }                             }
                            else { 
                                if (!pDataArray->keepprimer)    { currSeq.setAligned(currSeq.getUnaligned().substr(0, primerStart));   } 
                                else                { currSeq.setAligned(currSeq.getUnaligned().substr(0, primerEnd));     }
                            }
                        }
                    }
                }else if (pDataArray->ecolifile != "") {
                    //make sure the seqs are aligned
                    if (!pDataArray->fileAligned) { pDataArray->m->mothurOut("[ERROR]: seqs are not aligned. When using start and end your sequences must be aligned.\n"); pDataArray->m->setControl_pressed(true); break; }
                    else if (currSeq.getAligned().length() != pDataArray->length) {
                        pDataArray->m->mothurOut("[ERROR]: seqs are not the same length as ecoli seq. When using ecoli option your sequences must be aligned and the same length as the ecoli sequence.\n"); pDataArray->m->setControl_pressed(true); break; 
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
                    if (!pDataArray->fileAligned) { pDataArray->m->mothurOut("[ERROR]: seqs are not aligned. When using start and end your sequences must be aligned.\n"); pDataArray->m->setControl_pressed(true); break; }
                    else {
                        if (pDataArray->end != -1) {
                            if (pDataArray->end > currSeq.getAligned().length()) {  pDataArray->m->mothurOut("[ERROR]: end is longer than your sequence length, aborting.\n"); pDataArray->m->setControl_pressed(true); break; }
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
                
                if (commentString != "") {
                    string seqComment = currSeq.getComment();
                    currSeq.setComment("\t" + commentString + "\t" + seqComment);
                }
                
                if (totalDiffs > (pDataArray->pdiffs + pDataArray->rdiffs)) { trashCode += "t"; goodSeq = false; }
                
                //trimming removed all bases
                if (currSeq.getUnaligned() == "") { goodSeq = false; }
                
				if(goodSeq == 1)    {
                    currSeq.printSequence(goodFile);
                    if (locationsString != "") { locationsFile << locationsString; }
                    if (thisPStart != -1)   { locations.insert(thisPStart);  }
                }
				else {  
                    pDataArray->badSeqNames.insert(currSeq.getName()); 
                    currSeq.setName(currSeq.getName() + '|' + trashCode);
                    currSeq.printSequence(badFile); 
                }
			}
						
			//report progress
			if((i+1) % 100 == 0){	pDataArray->m->mothurOutJustToScreen("Processing sequence: " + toString(i+1)+"\n"); 	}
		}
		//report progress
		if((pDataArray->count) % 100 != 0){	pDataArray->m->mothurOutJustToScreen("Thread Processing sequence: " + toString(pDataArray->count)+"\n"); 		}
		
		goodFile.close();
		inFASTA.close();
        badFile.close();
        locationsFile.close();
        
        if (pDataArray->m->getDebug()) { pDataArray->m->mothurOut("[DEBUG]: fileAligned = " + toString(pDataArray->fileAligned) +'\n'); }
        
        if (pDataArray->fileAligned && !pDataArray->keepdots) { //print out smallest start value and largest end value
            if (locations.size() > 1) { pDataArray->adjustNeeded = true; }
            if (numFPrimers != 0)    {   set<int>::iterator it = locations.begin();  pDataArray->pstart = *it;  }
        }
        
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
