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
#include "oligos.h"
#include "fastqread.h"

struct pairFastqRead {
	FastqRead forward;
    FastqRead reverse;
    FastqRead findex;
    FastqRead rindex;
	
	pairFastqRead() {};
	pairFastqRead(FastqRead f, FastqRead r) : forward(f), reverse(r){};
    pairFastqRead(FastqRead f, FastqRead r, FastqRead fi, FastqRead ri) : forward(f), reverse(r), findex(fi), rindex(ri) {};
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
    struct linePair {
        unsigned long long start;
        unsigned long long end;
        linePair(unsigned long long i, unsigned long long j) : start(i), end(j) {}
        linePair() {}
    };
    
    char delim; 
    bool abort, allFiles, trimOverlap, createFileGroup, createOligosGroup, makeCount, noneOk, reorient;
    string outputDir, ffastqfile, rfastqfile, align, oligosfile, rfastafile, ffastafile, rqualfile, fqualfile, findexfile, rindexfile, file, format, inputDir;
    string outFastaFile, outQualFile, outScrapFastaFile, outScrapQualFile, outMisMatchFile, outputGroupFileName;
	float match, misMatch, gapOpen, gapExtend;
	int processors, longestBase, insert, tdiffs, bdiffs, pdiffs, ldiffs, sdiffs, deltaq, numBarcodes, numFPrimers, numLinkers, numSpacers, numRPrimers;
    vector<string> outputNames;
    Oligos* oligos;
    
	map<string, int> groupCounts; 
    map<string, string> groupMap;
    map<int, string> file2Group;
    
    unsigned long long processMultipleFileOption(map<string, int>&);
    unsigned long long processSingleFileOption(map<string, int>&);
    
    //main processing functions
    unsigned long long createProcesses(vector<string>, vector<string>, string, string, string, string, string, vector<vector<string> >, vector<vector<string> >, vector<linePair>, vector<linePair>, string);
    int driver(vector<string> files, vector<string> qualOrIndexFiles, string outputFasta, string outputScrapFasta, string outputQual, string outputScrapQual,  string outputMisMatches, vector<vector<string> > fastaFileNames, vector<vector<string> > qualFileNames, linePair, linePair, linePair, linePair, string);
    
    vector< vector<string> > readFileNames(string);
    bool getOligos(vector<vector<string> >&, vector<vector<string> >&, string, map<string, string>&);
    int setLines(vector<string>, vector<string>, vector<linePair>& fastaFilePos, vector<linePair>& qfileFilePos, char delim); //the delim let you know whether this is fasta and qual, or fastq and index. linePair entries will always be in sets of two. One for the forward and one for hte reverse.  (fastaFilePos[0] - ffasta, fastaFilePos[1] - rfasta) - processor1
};

/**************************************************************************************************/

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct contigsData {
    struct linePair {
        unsigned long long start;
        unsigned long long end;
        linePair(unsigned long long i, unsigned long long j) : start(i), end(j) {}
        linePair() {}
    };
    
    char delim;
    linePair linesInput, linesInputReverse, qlinesInput, qlinesInputReverse;
	string outputFasta, outputQual;
    string outputScrapFasta, outputScrapQual;
	string outputMisMatches;
	string align, group, oligosfile, format;
    vector<string> inputFiles, qualOrIndexFiles;
    vector<vector<string> > fastaFileNames, qualFileNames;
	MothurOut* m;
	float match, misMatch, gapOpen, gapExtend;
	int count, insert, threadID, pdiffs, bdiffs, tdiffs, deltaq;
    bool allFiles, createOligosGroup, createFileGroup, done, trimOverlap, reorient;
    map<string, int> groupCounts; 
    map<string, string> groupMap;
    
	
	contigsData(){}
	contigsData(string form, char d, string g, vector<string> f, vector<string> qif, string of, string osf, string oq, string osq, string om, string al, MothurOut* mout, float ma, float misMa, float gapO, float gapE, int thr, int delt, vector<vector<string> > ffn, vector<vector<string> > qfn,string olig, bool ro, int pdf, int bdf, int tdf, bool cg, bool cfg, bool all, bool to, linePair lff, linePair lrf, linePair qff, linePair qrf, int tid) {
        inputFiles = f;
        qualOrIndexFiles = qif;
		outputFasta = of;
        outputMisMatches = om;
        outputQual = oq;
        outputScrapQual = osq;
        m = mout;
		match = ma; 
		misMatch = misMa;
		gapOpen = gapO; 
		gapExtend = gapE; 
        insert = thr;
		align = al;
        group = g;
		count = 0;
        outputScrapFasta = osf;
        fastaFileNames = ffn;
        qualFileNames = qfn;
        oligosfile = olig;
        pdiffs = pdf;
        bdiffs = bdf;
        tdiffs = tdf;
        allFiles = all;
        trimOverlap = to;
        createOligosGroup = cg;
        createFileGroup = cfg;
		threadID = tid;
        deltaq = delt;
        reorient = ro;
        linesInput = lff;
        linesInputReverse = lrf;
        qlinesInput = qff;
        qlinesInputReverse = qrf;
        delim = d;
        format = form;
        done=false;
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
        
        pDataArray->count = 0;
        int num = 0;
        string thisfqualindexfile, thisrqualindexfile, thisffastafile, thisrfastafile;
        thisfqualindexfile = ""; thisrqualindexfile = "";
        thisffastafile = pDataArray->inputFiles[0]; thisrfastafile = pDataArray->inputFiles[1];
        if (pDataArray->qualOrIndexFiles.size() != 0) {
            thisfqualindexfile = pDataArray->qualOrIndexFiles[0];
            thisrqualindexfile = pDataArray->qualOrIndexFiles[1];
        }
        
        if (pDataArray->m->debug) {  pDataArray->m->mothurOut("[DEBUG]: ffasta = " + thisffastafile + ".\n[DEBUG]: rfasta = " + thisrfastafile + ".\n[DEBUG]: fqualindex = " + thisfqualindexfile + ".\n[DEBUG]: rqualindex = " + thisfqualindexfile + ".\n"); }
        
        ifstream inFFasta, inRFasta, inFQualIndex, inRQualIndex;
        ofstream outFasta, outMisMatch, outScrapFasta, outQual, outScrapQual;
        pDataArray->m->openInputFile(thisffastafile, inFFasta);
        pDataArray->m->openInputFile(thisrfastafile, inRFasta);
        
        inFFasta.seekg(pDataArray->linesInput.start);
        inRFasta.seekg(pDataArray->linesInputReverse.start);
        
        if (thisfqualindexfile != "") {
            if (thisfqualindexfile != "NONE") {
                pDataArray->m->openInputFile(thisfqualindexfile, inFQualIndex);
                inFQualIndex.seekg(pDataArray->qlinesInput.start);
            }
            else {  thisfqualindexfile = ""; }
            if (thisrqualindexfile != "NONE") {
                pDataArray->m->openInputFile(thisrqualindexfile, inRQualIndex);
                inRQualIndex.seekg(pDataArray->qlinesInputReverse.start);
            }
            else { thisrqualindexfile = ""; }
        }
        
        pDataArray->m->openOutputFile(pDataArray->outputFasta, outFasta);
        pDataArray->m->openOutputFile(pDataArray->outputScrapFasta, outScrapFasta);
        pDataArray->m->openOutputFile(pDataArray->outputMisMatches, outMisMatch);
        bool hasQuality = false;
        outMisMatch << "Name\tLength\tOverlap_Length\tOverlap_Start\tOverlap_End\tMisMatches\tNum_Ns\n";
        if (pDataArray->delim == '@') { //fastq files so make an output quality
            pDataArray->m->openOutputFile(pDataArray->outputQual, outQual);
            pDataArray->m->openOutputFile(pDataArray->outputScrapQual, outScrapQual);
            hasQuality = true;
        }else if ((pDataArray->delim == '>') && (pDataArray->qualOrIndexFiles.size() != 0)) { //fasta and qual files
            pDataArray->m->openOutputFile(pDataArray->outputQual, outQual);
            pDataArray->m->openOutputFile(pDataArray->outputScrapQual, outScrapQual);
            hasQuality = true;
        }
        
		if(pDataArray->allFiles){
			for (int i = 0; i < pDataArray->fastaFileNames.size(); i++) { //clears old file
				for (int j = 0; j < pDataArray->fastaFileNames[i].size(); j++) { //clears old file
					if (pDataArray->fastaFileNames[i][j] != "") {
						ofstream temp, temp2;
						pDataArray->m->openOutputFile(pDataArray->fastaFileNames[i][j], temp);			temp.close();
                        pDataArray->m->openOutputFile(pDataArray->qualFileNames[i][j], temp2);			temp2.close();
					}
				}
			}
		}
       
        Oligos oligos;
        if (pDataArray->oligosfile != "") { oligos.read(pDataArray->oligosfile, false);  }
        int numFPrimers = oligos.getPairedPrimers().size();
        int numBarcodes = oligos.getPairedBarcodes().size();

        
        TrimOligos trimOligos(pDataArray->pdiffs, pDataArray->bdiffs, 0, 0, oligos.getPairedPrimers(), oligos.getPairedBarcodes());
        TrimOligos* rtrimOligos = NULL;
        if (pDataArray->reorient) {
            rtrimOligos = new TrimOligos(pDataArray->pdiffs, pDataArray->bdiffs, 0, 0, oligos.getReorientedPairedPrimers(), oligos.getReorientedPairedBarcodes()); numBarcodes = oligos.getReorientedPairedBarcodes().size();
        }
        
        while ((!inFFasta.eof()) && (!inRFasta.eof())) {
            
            if (pDataArray->m->control_pressed) { break; }
            
            int success = 1;
            string trashCode = "";
            string commentString = "";
            int currentSeqsDiffs = 0;
            bool hasIndex = false;
            
            bool ignore; ignore = false;
            Sequence fSeq, rSeq;
            QualityScores* fQual = NULL; QualityScores* rQual = NULL;
            QualityScores* savedFQual = NULL; QualityScores* savedRQual = NULL;
            Sequence findexBarcode("findex", "NONE");  Sequence rindexBarcode("rindex", "NONE");
            if (pDataArray->delim == '@') { //fastq files
                bool tignore;
                FastqRead fread(inFFasta, tignore, pDataArray->format); pDataArray->m->gobble(inFFasta);
                FastqRead rread(inRFasta, ignore, pDataArray->format); pDataArray->m->gobble(inRFasta);
                if (tignore) { ignore=true; }
                fSeq.setName(fread.getName()); fSeq.setAligned(fread.getSeq());
                rSeq.setName(rread.getName()); rSeq.setAligned(rread.getSeq());
                fQual = new QualityScores(fread.getName(), fread.getScores());
                rQual = new QualityScores(rread.getName(), rread.getScores());
                if (thisfqualindexfile != "") { //forward index file
                    FastqRead firead(inFQualIndex, tignore, pDataArray->format); pDataArray->m->gobble(inFQualIndex);
                    if (tignore) { ignore=true; }
                    findexBarcode.setAligned(firead.getSeq());
                    if (firead.getName() != fread.getName()) { pDataArray->m->mothurOut("[WARNING]: name mismatch in forward index file. Ignoring, " + fread.getName() + ".\n"); ignore = true; }
                    hasIndex = true;
                }
                if (thisrqualindexfile != "") { //reverse index file
                    FastqRead riread(inRQualIndex, tignore, pDataArray->format); pDataArray->m->gobble(inRQualIndex);
                    if (tignore) { ignore=true; }
                    rindexBarcode.setAligned(riread.getSeq());
                    if (riread.getName() != fread.getName()) { pDataArray->m->mothurOut("[WARNING]: name mismatch in reverse index file. Ignoring, " + fread.getName() + ".\n"); ignore = true; }
                    hasIndex = true;
                }
                if (fread.getName() != rread.getName()) { pDataArray->m->mothurOut("[WARNING]: name mismatch in forward and reverse fastq file. Ignoring, " + fread.getName() + ".\n"); ignore = true; }
            }else { //reading fasta and maybe qual
                Sequence tfSeq(inFFasta); pDataArray->m->gobble(inFFasta);
                Sequence trSeq(inRFasta); pDataArray->m->gobble(inRFasta);
                fSeq.setName(tfSeq.getName()); fSeq.setAligned(tfSeq.getAligned());
                rSeq.setName(trSeq.getName()); rSeq.setAligned(trSeq.getAligned());
                if (thisfqualindexfile != "") {
                    fQual = new QualityScores(inFQualIndex); pDataArray->m->gobble(inFQualIndex);
                    rQual = new QualityScores(inRQualIndex); pDataArray->m->gobble(inRQualIndex);
                    savedFQual = new QualityScores(fQual->getName(), fQual->getQualityScores());
                    savedRQual = new QualityScores(rQual->getName(), rQual->getQualityScores());
                    if (fQual->getName() != tfSeq.getName()) { pDataArray->m->mothurOut("[WARNING]: name mismatch in forward quality file. Ignoring, " + tfSeq.getName() + ".\n"); ignore = true; }
                    if (rQual->getName() != trSeq.getName()) { pDataArray->m->mothurOut("[WARNING]: name mismatch in reverse quality file. Ignoring, " + trSeq.getName() + ".\n"); ignore = true; }
                }
                if (tfSeq.getName() != trSeq.getName()) { pDataArray->m->mothurOut("[WARNING]: name mismatch in forward and reverse fasta file. Ignoring, " + tfSeq.getName() + ".\n"); ignore = true; }
            }
            
            int barcodeIndex = 0;
            int primerIndex = 0;
            
            if (!ignore) {
                
                Sequence savedFSeq(fSeq.getName(), fSeq.getAligned());  Sequence savedRSeq(rSeq.getName(), rSeq.getAligned());
                Sequence savedFindex(findexBarcode.getName(), findexBarcode.getAligned()); Sequence savedRIndex(rindexBarcode.getName(), rindexBarcode.getAligned());
                
                if(numBarcodes != 0){
                    vector<int> results;
                    if (hasQuality) {
                        if (hasIndex) {
                            results = trimOligos.stripBarcode(findexBarcode, rindexBarcode, *fQual, *rQual, barcodeIndex);
                        }else {
                            results = trimOligos.stripBarcode(fSeq, rSeq, *fQual, *rQual, barcodeIndex);
                        }
                    }else {
                        results = trimOligos.stripBarcode(fSeq, rSeq, barcodeIndex);
                    }
                    success = results[0] + results[2];
                    commentString += "fbdiffs=" + toString(results[0]) + "(" + trimOligos.getCodeValue(results[1], pDataArray->bdiffs) + "), rbdiffs=" + toString(results[2]) + "(" + trimOligos.getCodeValue(results[3], pDataArray->bdiffs) + ") ";
                    if(success > pDataArray->bdiffs)		{	trashCode += 'b';	}
                    else{ currentSeqsDiffs += success;  }
                }
                
                if(numFPrimers != 0){
                    vector<int> results;
                    if (hasQuality) {
                        results = trimOligos.stripForward(fSeq, rSeq, *fQual, *rQual, primerIndex);
                    }else {
                        results = trimOligos.stripForward(fSeq, rSeq, primerIndex);
                    }
                    success = results[0] + results[2];
                    commentString += "fpdiffs=" + toString(results[0]) + "(" + trimOligos.getCodeValue(results[1], pDataArray->pdiffs) + "), rpdiffs=" + toString(results[2]) + "(" + trimOligos.getCodeValue(results[3], pDataArray->pdiffs) + ") ";
                    if(success > pDataArray->pdiffs)		{	trashCode += 'f';	}
                    else{ currentSeqsDiffs += success;  }
                }
                
                if (currentSeqsDiffs > pDataArray->tdiffs)	{	trashCode += 't';   }
                
                if (pDataArray->reorient && (trashCode != "")) { //if you failed and want to check the reverse
                    int thisSuccess = 0;
                    string thisTrashCode = "";
                    string thiscommentString = "";
                    int thisCurrentSeqsDiffs = 0;
                    
                    int thisBarcodeIndex = 0;
                    int thisPrimerIndex = 0;
                    
                    if(numBarcodes != 0){
                        vector<int> results;
                        if (hasQuality) {
                            if (hasIndex) {
                                results = rtrimOligos->stripBarcode(savedFindex, savedRIndex, *savedFQual, *savedRQual, thisBarcodeIndex);
                            }else {
                                results = rtrimOligos->stripBarcode(savedFSeq, savedRSeq, *savedFQual, *savedRQual, thisBarcodeIndex);
                            }
                        }else {
                            results = rtrimOligos->stripBarcode(savedFSeq, savedRSeq, thisBarcodeIndex);
                        }
                        thisSuccess = results[0] + results[2];
                        thiscommentString += "fbdiffs=" + toString(results[0]) + "(" + rtrimOligos->getCodeValue(results[1], pDataArray->bdiffs) + "), rbdiffs=" + toString(results[2]) + "(" + rtrimOligos->getCodeValue(results[3], pDataArray->bdiffs) + ") ";
                        if(thisSuccess > pDataArray->bdiffs)		{	thisTrashCode += 'b';	}
                        else{ thisCurrentSeqsDiffs += thisSuccess;  }
                    }
                    
                    if(numFPrimers != 0){
                        vector<int> results;
                        if (hasQuality) {
                            results = rtrimOligos->stripForward(savedFSeq, savedRSeq, *savedFQual, *savedRQual, thisPrimerIndex);
                        }else {
                            results = rtrimOligos->stripForward(savedFSeq, savedRSeq, thisPrimerIndex);
                        }
                        thisSuccess = results[0] + results[2];
                        thiscommentString += "fpdiffs=" + toString(results[0]) + "(" + rtrimOligos->getCodeValue(results[1], pDataArray->pdiffs) + "), rpdiffs=" + toString(results[2]) + "(" + rtrimOligos->getCodeValue(results[3], pDataArray->pdiffs) + ") ";
                        if(thisSuccess > pDataArray->pdiffs)		{	thisTrashCode += 'f';	}
                        else{ thisCurrentSeqsDiffs += thisSuccess;  }
                    }
                    
                    if (thisCurrentSeqsDiffs > pDataArray->tdiffs)	{	thisTrashCode += 't';   }
                    
                    if (thisTrashCode == "") {
                        trashCode = thisTrashCode;
                        success = thisSuccess;
                        currentSeqsDiffs = thisCurrentSeqsDiffs;
                        commentString = thiscommentString;
                        barcodeIndex = thisBarcodeIndex;
                        primerIndex = thisPrimerIndex;
                        savedFSeq.reverseComplement();
                        savedRSeq.reverseComplement();
                        fSeq.setAligned(savedFSeq.getAligned());
                        rSeq.setAligned(savedRSeq.getAligned());
                        if(hasQuality){
                            savedFQual->flipQScores(); savedRQual->flipQScores();
                            fQual->setScores(savedFQual->getScores()); rQual->setScores(savedRQual->getScores());
                        }
                    }else { trashCode += "(" + thisTrashCode + ")";  }
                }
                
                
                //flip the reverse reads
                rSeq.reverseComplement();
                if (hasQuality) { rQual->flipQScores(); }
                
                //pairwise align
                alignment->align(fSeq.getUnaligned(), rSeq.getUnaligned());
                map<int, int> ABaseMap = alignment->getSeqAAlnBaseMap();
                map<int, int> BBaseMap = alignment->getSeqBAlnBaseMap();
                fSeq.setAligned(alignment->getSeqAAln());
                rSeq.setAligned(alignment->getSeqBAln());
                int length = fSeq.getAligned().length();
                
                //traverse alignments merging into one contiguous seq
                string contig = "";
                int numMismatches = 0;
                string seq1 = fSeq.getAligned();
                string seq2 = rSeq.getAligned();
                vector<int> scores1, scores2, contigScores;
                if (hasQuality) {
                    scores1 = fQual->getQualityScores();
                    scores2 = rQual->getQualityScores();
                    delete fQual; delete rQual;  delete savedFQual; delete savedRQual;
                }
                
                // if (num < 5) {  cout << fSeq.getStartPos() << '\t' << fSeq.getEndPos() << '\t' << rSeq.getStartPos() << '\t' << rSeq.getEndPos() << endl; }
                int overlapStart = fSeq.getStartPos();
                int seq2Start = rSeq.getStartPos();
                
                //bigger of the 2 starting positions is the location of the overlapping start
                if (overlapStart < seq2Start) { //seq2 starts later so take from 0 to seq2Start from seq1
                    overlapStart = seq2Start;
                    for (int i = 0; i < overlapStart; i++) { contig += seq1[i];  }
                }else { //seq1 starts later so take from 0 to overlapStart from seq2
                    for (int i = 0; i < overlapStart; i++) {  contig += seq2[i]; }
                }
                
                int seq1End = fSeq.getEndPos();
                int seq2End = rSeq.getEndPos();
                int overlapEnd = seq1End;
                if (seq2End < overlapEnd) { overlapEnd = seq2End; }  //smallest end position is where overlapping ends
                
                int oStart = contig.length();
                //cout << fSeq.getAligned()  << endl; cout << rSeq.getAligned() << endl;
                for (int i = overlapStart; i < overlapEnd; i++) {
                    //cout << seq1[i] << ' ' << seq2[i] << ' ' << scores1[ABaseMap[i]] << ' ' << scores2[BBaseMap[i]] << endl;
                    if (seq1[i] == seq2[i]) { //match, add base and choose highest score
                        contig += seq1[i];
                    }else if (((seq1[i] == '.') || (seq1[i] == '-')) && ((seq2[i] != '-') && (seq2[i] != '.'))) { //seq1 is a gap and seq2 is a base, choose seq2, unless quality score for base is below insert. In that case eliminate base
                        if (hasQuality) {
                            if (scores2[BBaseMap[i]] <= pDataArray->insert) { } //
                            else { contig += seq2[i];  }
                        }else { contig += seq2[i]; } //with no quality info, then we keep it?
                    }else if (((seq2[i] == '.') || (seq2[i] == '-')) && ((seq1[i] != '-') && (seq1[i] != '.'))) { //seq2 is a gap and seq1 is a base, choose seq1, unless quality score for base is below insert. In that case eliminate base
                        if (hasQuality) {
                            if (scores1[ABaseMap[i]] <= pDataArray->insert) { } //
                            else { contig += seq1[i];  }
                        }else { contig += seq1[i]; } //with no quality info, then we keep it?
                    }else if (((seq1[i] != '-') && (seq1[i] != '.')) && ((seq2[i] != '-') && (seq2[i] != '.'))) { //both bases choose one with better quality
                        if (hasQuality) {
                            if (abs(scores1[ABaseMap[i]] - scores2[BBaseMap[i]]) >= pDataArray->deltaq) { //is the difference in qual scores >= deltaq, if yes choose base with higher score
                                char c = seq1[i];
                                if (scores1[ABaseMap[i]] < scores2[BBaseMap[i]]) { c = seq2[i]; }
                                contig += c;
                            }else { //if no, base becomes n
                                contig += 'N';
                            }
                            numMismatches++;
                        }else { numMismatches++; } //cant decide, so eliminate and mark as mismatch
                    }else { //should never get here
                        pDataArray->m->mothurOut("[ERROR]: case I didn't think of seq1 = " + toString(seq1[i]) + " and seq2 = " + toString(seq2[i]) + "\n");
                    }
                }
                int oend = contig.length();
                if (seq1End < seq2End) { //seq1 ends before seq2 so take from overlap to length from seq2
                    for (int i = overlapEnd; i < length; i++) { contig += seq2[i];  }
                }else { //seq2 ends before seq1 so take from overlap to length from seq1
                    for (int i = overlapEnd; i < length; i++) {  contig += seq1[i]; }
                }
                //cout << contig << endl;
                //exit(1);
                if (pDataArray->trimOverlap) { contig = contig.substr(overlapStart-1, oend-oStart);  if (contig.length() == 0) { trashCode += "l"; } }
                
                if(trashCode.length() == 0){
                    bool ignore = false;
                    
                    if (pDataArray->m->debug) { pDataArray->m->mothurOut(fSeq.getName()); }
                    
                    if (pDataArray->createOligosGroup) {
                        string thisGroup = oligos.getGroupName(barcodeIndex, primerIndex);
                        if (pDataArray->m->debug) { pDataArray->m->mothurOut(", group= " + thisGroup + "\n"); }
                        
                        int pos = thisGroup.find("ignore");
                        if (pos == string::npos) {
                            pDataArray->groupMap[fSeq.getName()] = thisGroup;
                            
                            map<string, int>::iterator it = pDataArray->groupCounts.find(thisGroup);
                            if (it == pDataArray->groupCounts.end()) {	pDataArray->groupCounts[thisGroup] = 1; }
                            else { pDataArray->groupCounts[it->first] ++; }
                        }else { ignore = true; }
                    }else if (pDataArray->createFileGroup) { //for 3 column file option
                        int pos = pDataArray->group.find("ignore");
                        if (pos == string::npos) {
                            pDataArray->groupMap[fSeq.getName()] = pDataArray->group;
                            
                            map<string, int>::iterator it = pDataArray->groupCounts.find(pDataArray->group);
                            if (it == pDataArray->groupCounts.end()) {	pDataArray->groupCounts[pDataArray->group] = 1; }
                            else { pDataArray->groupCounts[it->first] ++; }
                        }else { ignore = true; }
                    }
                    if (pDataArray->m->debug) { pDataArray->m->mothurOut("\n"); }
                    
                    if(!ignore){
                        //output
                        outFasta << ">" << fSeq.getName() << '\t' << commentString << endl << contig << endl;
                        outQual << ">" << fSeq.getName() << '\t' << commentString << endl;
                        for (int i = 0; i < contigScores.size(); i++) { outQual << contigScores[i] << " "; }  outQual << endl;
                        
                        int numNs = 0;
                        for (int i = 0; i < contig.length(); i++) { if (contig[i] == 'N') { numNs++; }  }
                        outMisMatch << fSeq.getName() << '\t' << contig.length() << '\t' << (oend-oStart) << '\t' << oStart << '\t' << oend << '\t' << numMismatches << '\t' << numNs << endl;
                        
                        if (pDataArray->allFiles) {
                            ofstream output;
                            pDataArray->m->openOutputFileAppend(pDataArray->fastaFileNames[barcodeIndex][primerIndex], output);
                            output << ">" << fSeq.getName() << '\t' << commentString << endl << contig << endl;
                            output.close();
                            
                            ofstream output2;
                            m->openOutputFileAppend(qualFileNames[barcodeIndex][primerIndex], output2);
                            output2 << ">" << fSeq.getName() << '\t' << commentString << endl;
                            for (int i = 0; i < contigScores.size(); i++) { output2 << contigScores[i] << " "; }  output2 << endl;
                            output2.close();

                        }
                    }
                }else {
                    //output
                    outScrapFasta << ">" << fSeq.getName() << " | " << trashCode << '\t' << commentString << endl << contig << endl;
                    outScrapQual << ">" << fSeq.getName() << '\t' << commentString << endl;
                    for (int i = 0; i < contigScores.size(); i++) { outScrapQual << contigScores[i] << " "; }  outScrapQual << endl;
                }
            }
            pDataArray->count++;
            
			//report progress
			if((pDataArray->count) % 1000 == 0){	pDataArray->m->mothurOut(toString(pDataArray->count)); pDataArray->m->mothurOutEndLine();		}
		}
        
		//report progress
		if((pDataArray->count) % 1000 != 0){	pDataArray->m->mothurOut(toString(pDataArray->count)); pDataArray->m->mothurOutEndLine();		}
        
        inFFasta.close();
        inRFasta.close();
        outFasta.close();
        outScrapFasta.close();
        outMisMatch.close();
        if (pDataArray->delim == '@') {
            if (thisfqualindexfile != "") { inFQualIndex.close(); }
            if (thisrqualindexfile != "") { inRQualIndex.close(); }
            outQual.close();
            outScrapQual.close();
        }else{
            if (hasQuality) {
                inFQualIndex.close();
                inRQualIndex.close();
                outQual.close();
                outScrapQual.close();
            }
        }
        delete alignment;
        if (pDataArray->reorient) { delete rtrimOligos; }
        
        pDataArray->done = true;
        if (pDataArray->m->control_pressed) {  pDataArray->m->mothurRemove(pDataArray->outputFasta);  pDataArray->m->mothurRemove(pDataArray->outputMisMatches);  pDataArray->m->mothurRemove(pDataArray->outputScrapFasta); }
        
        return 0;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "AlignCommand", "MyContigsThreadFunction");
		exit(1);
	}
} 
#endif


#endif
