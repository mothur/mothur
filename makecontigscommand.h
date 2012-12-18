#ifndef Mothur_makecontigscommand_h
#define Mothur_makecontigscommand_h

//
//  makecontigscommand.h
//  Mothur
//
//  Created by Sarah Westcott on 5/15/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "command.hpp"
#include "sequence.hpp"
#include "qualityscores.h"
#include "alignment.hpp"
#include "gotohoverlap.hpp"
#include "needlemanoverlap.hpp"
#include "blastalign.hpp"
#include "noalign.hpp"
#include "trimoligos.h"

struct fastqRead {
	vector<int> scores;
	string name;
	string sequence;
	
	fastqRead() { name = ""; sequence = ""; scores.clear(); };
	fastqRead(string n, string s, vector<int> sc) : name(n), sequence(s), scores(sc) {};
	~fastqRead() {};
};

struct pairFastqRead {
	fastqRead forward;
    fastqRead reverse;
	
	pairFastqRead() {};
	pairFastqRead(fastqRead f, fastqRead r) : forward(f), reverse(r){};
	~pairFastqRead() {};
};
/**************************************************************************************************/

class MakeContigsCommand : public Command {
public:
    MakeContigsCommand(string);
    MakeContigsCommand();
    ~MakeContigsCommand(){}
    
    vector<string> setParameters();
    string getCommandName()			{ return "make.contigs";			}
    string getCommandCategory()		{ return "Sequence Processing";		} 
    //commmand category choices: Sequence Processing, OTU-Based Approaches, Hypothesis Testing, Phylotype Analysis, General, Clustering and Hidden
    
	string getHelpString();	
    string getOutputPattern(string);	
    string getCitation() { return "http://www.mothur.org/wiki/Make.contigs"; }
    string getDescription()		{ return "description"; }
    
    int execute(); 
    void help() { m->mothurOut(getHelpString()); }	
    
private:
    bool abort, allFiles, createGroup;
    string outputDir, ffastqfile, rfastqfile, align, oligosfile, rfastafile, ffastafile, rqualfile, fqualfile, file;
	float match, misMatch, gapOpen, gapExtend;
	int processors, longestBase, threshold, tdiffs, bdiffs, pdiffs, ldiffs, sdiffs;
    vector<string> outputNames;
    
    map<int, oligosPair> barcodes;
	map<int, oligosPair> primers;
    vector<string>  linker;
    vector<string>  spacer;
	vector<string> primerNameVector;	
	vector<string> barcodeNameVector;	
    
	map<string, int> groupCounts; 
    map<string, string> groupMap;
    //map<string, int> combos;
	//map<string, int> groupToIndex;
    //vector<string> groupVector;
    
    fastqRead readFastq(ifstream&, bool&);
    vector< vector<string> > readFastqFiles(unsigned long int&);
    bool checkReads(fastqRead&, fastqRead&);
    int createProcesses(vector< vector<string> >, string, string, string, string, string, vector<vector<string> >, vector<vector<string> >);
    int driver(vector<string>, string, string, string, string, string, vector<vector<string> >, vector<vector<string> >);
    bool getOligos(vector<vector<string> >&, vector<vector<string> >&);
    string reverseOligo(string);
    vector<pairFastqRead> getReads(bool ignoref, bool ignorer, fastqRead forward, fastqRead reverse, map<string, fastqRead>& uniques);
};

/**************************************************************************************************/

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct contigsData {
	string outputFasta; 
	string outputQual; 
    string outputScrapFasta; 
	string outputScrapQual;
	string outputMisMatches;
	string align;
    vector<string> files;
    vector<vector<string> > fastaFileNames;
    vector<vector<string> > qualFileNames;
	MothurOut* m;
	float match, misMatch, gapOpen, gapExtend;
	int count, threshold, threadID, pdiffs, bdiffs, tdiffs;
    bool allFiles, createGroup;
    map<string, int> groupCounts; 
    map<string, string> groupMap;
    vector<string> primerNameVector;	
	vector<string> barcodeNameVector;
    map<int, oligosPair> barcodes;
	map<int, oligosPair> primers;
	
	contigsData(){}
	contigsData(vector<string> f, string of, string oq, string osf, string osq, string om, string al, MothurOut* mout, float ma, float misMa, float gapO, float gapE, int thr, map<int, oligosPair> br, map<int, oligosPair> pr, vector<vector<string> > ffn, vector<vector<string> > qfn, vector<string>bnv, vector<string> pnv, int pdf, int bdf, int tdf, bool cg, bool all, int tid) {
        files = f;
		outputFasta = of;
        outputQual = oq;
        outputMisMatches = om;
        m = mout;
		match = ma; 
		misMatch = misMa;
		gapOpen = gapO; 
		gapExtend = gapE; 
        threshold = thr;
		align = al;
		count = 0;
        outputScrapFasta = osf;
        outputScrapQual = osq;
        fastaFileNames = ffn;
        qualFileNames = qfn;
        barcodes = br;
        primers = pr;
        barcodeNameVector = bnv;
        primerNameVector = pnv;
        pdiffs = pdf;
        bdiffs = bdf;
        tdiffs = tdf;
        allFiles = all;
        createGroup = cg;
		threadID = tid;
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyContigsThreadFunction(LPVOID lpParam){ 
	contigsData* pDataArray;
	pDataArray = (contigsData*)lpParam;
	
	try {
        int longestBase = 1000;
        Alignment* alignment;
        if(pDataArray->align == "gotoh")			{	alignment = new GotohOverlap(pDataArray->gapOpen, pDataArray->gapExtend, pDataArray->match, pDataArray->misMatch, longestBase);			}
		else if(pDataArray->align == "needleman")	{	alignment = new NeedlemanOverlap(pDataArray->gapOpen, pDataArray->match, pDataArray->misMatch, longestBase);				}
        
        int num = 0;
        string thisffastafile = pDataArray->files[0];
        string thisfqualfile = pDataArray->files[1];
        string thisrfastafile = pDataArray->files[2];
        string thisrqualfile = pDataArray->files[3];
        
        if (pDataArray->m->debug) {  pDataArray->m->mothurOut("[DEBUG]: ffasta = " + thisffastafile + ".\n[DEBUG]: fqual = " + thisfqualfile + ".\n[DEBUG]: rfasta = " + thisrfastafile + ".\n[DEBUG]: rqual = " + thisrqualfile + ".\n"); }
        
		if(pDataArray->allFiles){
			for (int i = 0; i < pDataArray->fastaFileNames.size(); i++) { //clears old file
				for (int j = 0; j < pDataArray->fastaFileNames[i].size(); j++) { //clears old file
					if (pDataArray->fastaFileNames[i][j] != "") {
						ofstream temp;
						pDataArray->m->openOutputFile(pDataArray->fastaFileNames[i][j], temp);			temp.close();
                        pDataArray->m->openOutputFile(pDataArray->qualFileNames[i][j], temp);			temp.close();
					}
				}
			}
		}
        
        ifstream inFFasta, inRFasta, inFQual, inRQual;
        pDataArray->m->openInputFile(thisffastafile, inFFasta);
        pDataArray->m->openInputFile(thisfqualfile, inFQual);
        pDataArray->m->openInputFile(thisrfastafile, inRFasta);
        pDataArray->m->openInputFile(thisrqualfile, inRQual);
        
        ofstream outFasta, outQual, outMisMatch, outScrapFasta, outScrapQual;
        pDataArray->m->openOutputFile(pDataArray->outputFasta, outFasta);
        pDataArray->m->openOutputFile(pDataArray->outputQual, outQual);
        pDataArray->m->openOutputFile(pDataArray->outputMisMatches, outMisMatch);
        pDataArray->m->openOutputFile(pDataArray->outputScrapFasta, outScrapFasta);
        pDataArray->m->openOutputFile(pDataArray->outputScrapQual, outScrapQual);
        outMisMatch << "Name\tLength\tMisMatches\n";
        
        TrimOligos trimOligos(pDataArray->pdiffs, pDataArray->bdiffs, 0, 0, pDataArray->primers, pDataArray->barcodes);
        
        while ((!inFQual.eof()) && (!inFFasta.eof()) && (!inRFasta.eof()) && (!inRQual.eof())) {
            
            if (pDataArray->m->control_pressed) { break; }
            
            int success = 1;
            string trashCode = "";
            int currentSeqsDiffs = 0;
            
            //read seqs and quality info
            Sequence fSeq(inFFasta); pDataArray->m->gobble(inFFasta);
            Sequence rSeq(inRFasta); pDataArray->m->gobble(inRFasta);
            QualityScores fQual(inFQual); pDataArray->m->gobble(inFQual);
            QualityScores rQual(inRQual); pDataArray->m->gobble(inRQual);
            
            int barcodeIndex = 0;
            int primerIndex = 0;
            
            if(pDataArray->barcodes.size() != 0){
                success = trimOligos.stripBarcode(fSeq, rSeq, fQual, rQual, barcodeIndex);
                if(success > pDataArray->bdiffs)		{	trashCode += 'b';	}
                else{ currentSeqsDiffs += success;  }
            }
            
            if(pDataArray->primers.size() != 0){
                success = trimOligos.stripForward(fSeq, rSeq, fQual, rQual, primerIndex);
                if(success > pDataArray->pdiffs)		{	trashCode += 'f';	}
                else{ currentSeqsDiffs += success;  }
            }
            
            if (currentSeqsDiffs > pDataArray->tdiffs)	{	trashCode += 't';   }
            
            //flip the reverse reads
            rSeq.reverseComplement();
            rQual.flipQScores();
           
            //pairwise align
            alignment->align(fSeq.getUnaligned(), rSeq.getUnaligned());
            map<int, int> ABaseMap = alignment->getSeqAAlnBaseMap();
            map<int, int> BBaseMap = alignment->getSeqBAlnBaseMap();
            fSeq.setAligned(alignment->getSeqAAln());
            rSeq.setAligned(alignment->getSeqBAln());
            int length = fSeq.getAligned().length();
            
            //traverse alignments merging into one contiguous seq
            string contig = "";
            vector<int> contigScores; 
            int numMismatches = 0;
            string seq1 = fSeq.getAligned();
            string seq2 = rSeq.getAligned();
            vector<int> scores1 = fQual.getQualityScores();
            vector<int> scores2 = rQual.getQualityScores();
            
            int overlapStart = fSeq.getStartPos();
            int seq2Start = rSeq.getStartPos();
            //bigger of the 2 starting positions is the location of the overlapping start
            if (overlapStart < seq2Start) { //seq2 starts later so take from 0 to seq2Start from seq1
                overlapStart = seq2Start; 
                for (int i = 0; i < overlapStart; i++) {
                    contig += seq1[i];
                    contigScores.push_back(scores1[ABaseMap[i]]);
                }
            }else { //seq1 starts later so take from 0 to overlapStart from seq2
                for (int i = 0; i < overlapStart; i++) {
                    contig += seq2[i];
                    contigScores.push_back(scores2[BBaseMap[i]]);
                }
            }
            
            int seq1End = fSeq.getEndPos();
            int seq2End = rSeq.getEndPos();
            int overlapEnd = seq1End;
            if (seq2End < overlapEnd) { overlapEnd = seq2End; }  //smallest end position is where overlapping ends
            
            for (int i = overlapStart; i < overlapEnd; i++) {
                if (seq1[i] == seq2[i]) { //match, add base and choose highest score
                    contig += seq1[i];
                    contigScores.push_back(scores1[ABaseMap[i]]);
                    if (scores1[ABaseMap[i]] < scores2[BBaseMap[i]]) { contigScores[contigScores.size()-1] = scores2[BBaseMap[i]]; }
                }else if (((seq1[i] == '.') || (seq1[i] == '-')) && ((seq2[i] != '-') && (seq2[i] != '.'))) { //seq1 is a gap and seq2 is a base, choose seq2, unless quality score for base is below threshold. In that case eliminate base
                    if (scores2[BBaseMap[i]] < pDataArray->threshold) { } //
                    else {
                        contig += seq2[i];
                        contigScores.push_back(scores2[BBaseMap[i]]);
                    }
                }else if (((seq2[i] == '.') || (seq2[i] == '-')) && ((seq1[i] != '-') && (seq1[i] != '.'))) { //seq2 is a gap and seq1 is a base, choose seq1, unless quality score for base is below threshold. In that case eliminate base
                    if (scores1[ABaseMap[i]] < pDataArray->threshold) { } //
                    else {
                        contig += seq1[i];
                        contigScores.push_back(scores1[ABaseMap[i]]);
                    }
                }else if (((seq1[i] != '-') && (seq1[i] != '.')) && ((seq2[i] != '-') && (seq2[i] != '.'))) { //both bases choose one with better quality
                    char c = seq1[i];
                    contigScores.push_back(scores1[ABaseMap[i]]);
                    if (scores1[ABaseMap[i]] < scores2[BBaseMap[i]]) { contigScores[contigScores.size()-1] = scores2[BBaseMap[i]]; c = seq2[i]; }
                    contig += c;
                    numMismatches++;
                }else { //should never get here
                    pDataArray->m->mothurOut("[ERROR]: case I didn't think of seq1 = " + toString(seq1[i]) + " and seq2 = " + toString(seq2[i]) + "\n");
                }
            }
            
            if (seq1End < seq2End) { //seq1 ends before seq2 so take from overlap to length from seq2
                for (int i = overlapEnd; i < length; i++) {
                    contig += seq2[i];
                    contigScores.push_back(scores2[BBaseMap[i]]);
                }
            }else { //seq2 ends before seq1 so take from overlap to length from seq1
                for (int i = overlapEnd; i < length; i++) {
                    contig += seq1[i];
                    contigScores.push_back(scores1[ABaseMap[i]]);
                }
                
            }

            if(trashCode.length() == 0){
                if (pDataArray->createGroup) {
                    if(pDataArray->barcodes.size() != 0){
                        string thisGroup = pDataArray->barcodeNameVector[barcodeIndex];
                        if (pDataArray->primers.size() != 0) { 
                            if (pDataArray->primerNameVector[primerIndex] != "") { 
                                if(thisGroup != "") {
                                    thisGroup += "." + pDataArray->primerNameVector[primerIndex]; 
                                }else {
                                    thisGroup = pDataArray->primerNameVector[primerIndex]; 
                                }
                            } 
                        }
                        
                        if (pDataArray->m->debug) { pDataArray->m->mothurOut(", group= " + thisGroup + "\n"); }
                        
                        pDataArray->groupMap[fSeq.getName()] = thisGroup; 
                        
                        map<string, int>::iterator it = pDataArray->groupCounts.find(thisGroup);
                        if (it == pDataArray->groupCounts.end()) {	pDataArray->groupCounts[thisGroup] = 1; }
                        else { pDataArray->groupCounts[it->first] ++; }
                        
                    }
                }
                
                if(pDataArray->allFiles){
                    ofstream output;
                    pDataArray->m->openOutputFileAppend(pDataArray->fastaFileNames[barcodeIndex][primerIndex], output);
                    output << ">" << fSeq.getName() << endl << contig << endl;
                    output.close();
                    
                    pDataArray->m->openOutputFileAppend(pDataArray->qualFileNames[barcodeIndex][primerIndex], output);
                    output << ">" << fSeq.getName() << endl;
                    for (int i = 0; i < contigScores.size(); i++) { output << contigScores[i] << ' '; }
                    output << endl;
                    output.close();							
                }
                
                //output
                outFasta << ">" << fSeq.getName() << endl << contig << endl;
                outQual << ">" << fSeq.getName() << endl;
                for (int i = 0; i < contigScores.size(); i++) { outQual << contigScores[i] << ' '; }
                outQual << endl;
                outMisMatch << fSeq.getName() << '\t' << contig.length() << '\t' << numMismatches << endl;
            }else {
                //output
                outScrapFasta << ">" << fSeq.getName() << " | " << trashCode << endl << contig << endl;
                outScrapQual << ">" << fSeq.getName() << " | " << trashCode << endl;
                for (int i = 0; i < contigScores.size(); i++) { outScrapQual << contigScores[i] << ' '; }
                outScrapQual << endl;
            }
            num++;
            
			//report progress
			if((num) % 1000 == 0){	pDataArray->m->mothurOut(toString(num)); pDataArray->m->mothurOutEndLine();		}
		}
        
		//report progress
		if((num) % 1000 != 0){	pDataArray->m->mothurOut(toString(num)); pDataArray->m->mothurOutEndLine();		}
        
        inFFasta.close();
        inFQual.close();
        inRFasta.close();
        inRQual.close();
        outFasta.close();
        outQual.close();
        outMisMatch.close();
        outScrapFasta.close();
        outScrapQual.close();
        delete alignment;
        
        if (pDataArray->m->control_pressed) { pDataArray->m->mothurRemove(pDataArray->outputQual); pDataArray->m->mothurRemove(pDataArray->outputFasta);  pDataArray->m->mothurRemove(pDataArray->outputMisMatches); pDataArray->m->mothurRemove(pDataArray->outputScrapQual); pDataArray->m->mothurRemove(pDataArray->outputScrapFasta);}
        
        return 0;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "AlignCommand", "MyContigsThreadFunction");
		exit(1);
	}
} 
#endif


#endif
