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
    fastqRead findex;
    fastqRead rindex;
	
	pairFastqRead() {};
	pairFastqRead(fastqRead f, fastqRead r) : forward(f), reverse(r){};
    pairFastqRead(fastqRead f, fastqRead r, fastqRead fi, fastqRead ri) : forward(f), reverse(r), findex(fi), rindex(ri) {};
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
    bool abort, allFiles, trimOverlap, createFileGroup, createOligosGroup, makeCount, noneOk, reorient;
    string outputDir, ffastqfile, rfastqfile, align, oligosfile, rfastafile, ffastafile, rqualfile, fqualfile, findexfile, rindexfile, file, format, inputDir;
	float match, misMatch, gapOpen, gapExtend;
	int processors, longestBase, insert, tdiffs, bdiffs, pdiffs, ldiffs, sdiffs, deltaq, numBarcodes, numFPrimers, numLinkers, numSpacers, numRPrimers;
    vector<string> outputNames;
    Oligos* oligos;
	vector<char> convertTable;
    
	map<string, int> groupCounts; 
    map<string, string> groupMap;
    map<int, string> file2Group;
    
    vector<int> convertQual(string);
    fastqRead readFastq(ifstream&, bool&);
    vector< vector< vector<string> > > preProcessData(unsigned long int&);
    vector< vector<string> > readFileNames(string);
    vector< vector<string> > readFastqFiles(unsigned long int&, string, string, string, string);
    vector< vector<string> > readFastaFiles(unsigned long int&, string, string);
    //bool checkReads(fastqRead&, fastqRead&, string, string);
    int createProcesses(vector< vector<string> >, string, string, string, vector<vector<string> >, int);
    int driver(vector<string>, string, string, string, vector<vector<string> >, int, string);
    bool getOligos(vector<vector<string> >&, string, map<string, string>&);
    vector<pairFastqRead> getReads(bool ignoref, bool ignorer, fastqRead forward, fastqRead reverse, map<string, fastqRead>& uniques, bool);
    vector<pairFastqRead> mergeReads(vector<pairFastqRead> frReads, vector<pairFastqRead> friReads, map<string, pairFastqRead>& pairUniques);
};

/**************************************************************************************************/

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct contigsData {
	string outputFasta; 
    string outputScrapFasta; 
	string outputMisMatches;
	string align, group, oligosfile;
    vector<string> files;
    vector<vector<string> > fastaFileNames;
	MothurOut* m;
	float match, misMatch, gapOpen, gapExtend;
	int count, insert, threadID, pdiffs, bdiffs, tdiffs, deltaq;
    bool allFiles, createOligosGroup, createFileGroup, done, trimOverlap, reorient;
    map<string, int> groupCounts; 
    map<string, string> groupMap;
    
	
	contigsData(){}
	contigsData(string g, vector<string> f, string of, string osf, string om, string al, MothurOut* mout, float ma, float misMa, float gapO, float gapE, int thr, int delt, vector<vector<string> > ffn, string olig, bool ro, int pdf, int bdf, int tdf, bool cg, bool cfg, bool all, bool to, int tid) {
        files = f;
		outputFasta = of;
        outputMisMatches = om;
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
        string thisffastafile = pDataArray->files[0];
        string thisfqualfile = pDataArray->files[1];
        string thisrfastafile = pDataArray->files[2];
        string thisrqualfile = pDataArray->files[3];
        string thisfindexfile = pDataArray->files[4];
        string thisrindexfile = pDataArray->files[5];
        
        if (pDataArray->m->debug) {  pDataArray->m->mothurOut("[DEBUG]: ffasta = " + thisffastafile + ".\n[DEBUG]: fqual = " + thisfqualfile + ".\n[DEBUG]: rfasta = " + thisrfastafile + ".\n[DEBUG]: rqual = " + thisrqualfile + ".\n[DEBUG]: findex = " + thisfindexfile + ".\n[DEBUG]: rindex = " + thisrindexfile + ".\n"); }
        
		if(pDataArray->allFiles){
			for (int i = 0; i < pDataArray->fastaFileNames.size(); i++) { //clears old file
				for (int j = 0; j < pDataArray->fastaFileNames[i].size(); j++) { //clears old file
					if (pDataArray->fastaFileNames[i][j] != "") {
						ofstream temp;
						pDataArray->m->openOutputFile(pDataArray->fastaFileNames[i][j], temp);			temp.close();
					}
				}
			}
		}
        
        ifstream inFFasta, inRFasta, inFQual, inRQual, inFIndex, inRIndex;
        ofstream outFasta, outMisMatch, outScrapFasta;
        pDataArray->m->openInputFile(thisffastafile, inFFasta);
        pDataArray->m->openInputFile(thisrfastafile, inRFasta);
        if (thisfqualfile != "") {
            pDataArray->m->openInputFile(thisfqualfile, inFQual);
            pDataArray->m->openInputFile(thisrqualfile, inRQual);
        }
        
        if (thisfindexfile != "") { pDataArray->m->openInputFile(thisfindexfile, inFIndex);  }
        if (thisrindexfile != "") { pDataArray->m->openInputFile(thisrindexfile, inRIndex);  }
        
        pDataArray->m->openOutputFile(pDataArray->outputFasta, outFasta);
        pDataArray->m->openOutputFile(pDataArray->outputMisMatches, outMisMatch);
        pDataArray->m->openOutputFile(pDataArray->outputScrapFasta, outScrapFasta);
        
        outMisMatch << "Name\tLength\tOverlap_Length\tOverlap_Start\tOverlap_End\tMisMatches\tNum_Ns\n";  
        
        Oligos oligos;
        if (pDataArray->oligosfile != "") { oligos.read(pDataArray->oligosfile);  }
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
            int currentSeqsDiffs = 0;
            
            //read seqs and quality info
            Sequence fSeq(inFFasta); pDataArray->m->gobble(inFFasta);
            Sequence rSeq(inRFasta); pDataArray->m->gobble(inRFasta);
            QualityScores* fQual = NULL; QualityScores* rQual = NULL;
            if (thisfqualfile != "") {
                fQual = new QualityScores(inFQual); pDataArray->m->gobble(inFQual);
                rQual = new QualityScores(inRQual); pDataArray->m->gobble(inRQual);
            }
            
            Sequence findexBarcode("findex", "NONE");  Sequence rindexBarcode("rindex", "NONE");
            if (thisfindexfile != "") {
                Sequence temp(inFIndex); pDataArray->m->gobble(inFIndex);
                findexBarcode.setAligned(temp.getAligned());
            }
            
            if (thisrindexfile != "") {
                Sequence temp(inRIndex); pDataArray->m->gobble(inRIndex);
                rindexBarcode.setAligned(temp.getAligned());
            }

            int barcodeIndex = 0;
            int primerIndex = 0;
            Sequence savedFSeq(fSeq.getName(), fSeq.getAligned());  Sequence savedRSeq(rSeq.getName(), rSeq.getAligned());
            Sequence savedFindex(findexBarcode.getName(), findexBarcode.getAligned()); Sequence savedRIndex(rindexBarcode.getName(), rindexBarcode.getAligned());
            QualityScores* savedFQual = NULL; QualityScores* savedRQual = NULL;
            if (thisfqualfile != "") {
                savedFQual = new QualityScores(fQual->getName(), fQual->getQualityScores());
                savedRQual = new QualityScores(rQual->getName(), rQual->getQualityScores());
            }
            
            if(numBarcodes != 0){
                if (thisfqualfile != "") {
                    if ((thisfindexfile != "") || (thisrindexfile != "")) {
                        success = trimOligos.stripBarcode(findexBarcode, rindexBarcode, *fQual, *rQual, barcodeIndex);
                    }else {
                        success = trimOligos.stripBarcode(fSeq, rSeq, *fQual, *rQual, barcodeIndex);
                    }
                }else {
                    success = trimOligos.stripBarcode(fSeq, rSeq, barcodeIndex);
                }
                if(success > pDataArray->bdiffs)		{	trashCode += 'b';	}
                else{ currentSeqsDiffs += success;  }
            }
            
            if(numFPrimers != 0){
                if (thisfqualfile != "") {
                    success = trimOligos.stripForward(fSeq, rSeq, *fQual, *rQual, primerIndex);
                }else {
                    success = trimOligos.stripForward(fSeq, rSeq, primerIndex);
                }
                if(success > pDataArray->pdiffs)		{	trashCode += 'f';	}
                else{ currentSeqsDiffs += success;  }
            }
            
            if (currentSeqsDiffs > pDataArray->tdiffs)	{	trashCode += 't';   }
            
            if (pDataArray->reorient && (trashCode != "")) { //if you failed and want to check the reverse
                int thisSuccess = 0;
                string thisTrashCode = "";
                int thisCurrentSeqsDiffs = 0;
                
                int thisBarcodeIndex = 0;
                int thisPrimerIndex = 0;
                
                if(numBarcodes != 0){
                    if (thisfqualfile != "") {
                        if ((thisfindexfile != "") || (thisrindexfile != "")) {
                            thisSuccess = rtrimOligos->stripBarcode(savedFindex, savedRIndex, *savedFQual, *savedRQual, thisBarcodeIndex);
                        }else {
                            thisSuccess = rtrimOligos->stripBarcode(savedFSeq, savedRSeq, *savedFQual, *savedRQual, thisBarcodeIndex);
                        }
                    }else {
                        thisSuccess = rtrimOligos->stripBarcode(savedFSeq, savedRSeq, thisBarcodeIndex);
                    }
                    if(thisSuccess > pDataArray->bdiffs)		{	thisTrashCode += 'b';	}
                    else{ thisCurrentSeqsDiffs += thisSuccess;  }
                }
                
                if(numFPrimers != 0){
                    if (thisfqualfile != "") {
                        thisSuccess = rtrimOligos->stripForward(savedFSeq, savedRSeq, *savedFQual, *savedRQual, thisPrimerIndex);
                    }else {
                        thisSuccess = rtrimOligos->stripForward(savedFSeq, savedRSeq, thisPrimerIndex);
                    }
                    if(thisSuccess > pDataArray->pdiffs)		{	thisTrashCode += 'f';	}
                    else{ thisCurrentSeqsDiffs += thisSuccess;  }
                }
                
                if (thisCurrentSeqsDiffs > pDataArray->tdiffs)	{	thisTrashCode += 't';   }
                
                if (thisTrashCode == "") {
                    trashCode = thisTrashCode;
                    success = thisSuccess;
                    currentSeqsDiffs = thisCurrentSeqsDiffs;
                    barcodeIndex = thisBarcodeIndex;
                    primerIndex = thisPrimerIndex;
                    savedFSeq.reverseComplement();
                    savedRSeq.reverseComplement();
                    fSeq.setAligned(savedFSeq.getAligned());
                    rSeq.setAligned(savedRSeq.getAligned());
                    if(thisfqualfile != ""){
                        savedFQual->flipQScores(); savedRQual->flipQScores();
                        fQual->setScores(savedFQual->getScores()); rQual->setScores(savedRQual->getScores());
                    }
                }else { trashCode += "(" + thisTrashCode + ")";  }
            }
            
            //flip the reverse reads
            rSeq.reverseComplement();
            if (thisfqualfile != "") { rQual->flipQScores(); }
           
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
            vector<int> scores1, scores2;
            if (thisfqualfile != "") {
                scores1 = fQual->getQualityScores();
                scores2 = rQual->getQualityScores();
                delete fQual; delete rQual; delete savedFQual; delete savedRQual;
            }
            
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
            for (int i = overlapStart; i < overlapEnd; i++) {
                if (seq1[i] == seq2[i]) { //match, add base and choose highest score
                    contig += seq1[i];
                }else if (((seq1[i] == '.') || (seq1[i] == '-')) && ((seq2[i] != '-') && (seq2[i] != '.'))) { //seq1 is a gap and seq2 is a base, choose seq2, unless quality score for base is below insert. In that case eliminate base
                    if (thisfqualfile != "") {
                        if (scores2[BBaseMap[i]] < pDataArray->insert) { } //
                        else { contig += seq2[i];  }
                    }else { contig += seq2[i]; } //with no quality info, then we keep it?
                }else if (((seq2[i] == '.') || (seq2[i] == '-')) && ((seq1[i] != '-') && (seq1[i] != '.'))) { //seq2 is a gap and seq1 is a base, choose seq1, unless quality score for base is below insert. In that case eliminate base
                    if (thisfqualfile != "") {
                        if (scores1[ABaseMap[i]] < pDataArray->insert) { } //
                        else { contig += seq1[i];  }
                    }else { contig += seq1[i]; } //with no quality info, then we keep it?
                }else if (((seq1[i] != '-') && (seq1[i] != '.')) && ((seq2[i] != '-') && (seq2[i] != '.'))) { //both bases choose one with better quality
                    if (thisfqualfile != "") {
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

            if (pDataArray->trimOverlap) { contig = contig.substr(overlapStart-1, oend-oStart); if (contig.length() == 0) { trashCode += "l"; } }
            
            if(trashCode.length() == 0){
                bool ignore = false;
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
                }else if (pDataArray->createFileGroup) {
                    int pos = pDataArray->group.find("ignore");
                    if (pos == string::npos) {
                        pDataArray->groupMap[fSeq.getName()] = pDataArray->group;
                        
                        map<string, int>::iterator it = pDataArray->groupCounts.find(pDataArray->group);
                        if (it == pDataArray->groupCounts.end()) {	pDataArray->groupCounts[pDataArray->group] = 1; }
                        else { pDataArray->groupCounts[it->first]++; }
                    }else { ignore = true; }
                }

                
                if(pDataArray->allFiles && !ignore){
                    ofstream output;
                    pDataArray->m->openOutputFileAppend(pDataArray->fastaFileNames[barcodeIndex][primerIndex], output);
                    output << ">" << fSeq.getName() << endl << contig << endl;
                    output.close();
                }
                
                //output
                outFasta << ">" << fSeq.getName() << endl << contig << endl;
                int numNs = 0;
                for (int i = 0; i < contig.length(); i++) { if (contig[i] == 'N') { numNs++; }  }
                outMisMatch << fSeq.getName() << '\t' << contig.length() << '\t' << (oend-oStart) << '\t' << oStart << '\t' << oend << '\t' << numMismatches << '\t' << numNs << endl;
            }else {
                //output
                outScrapFasta << ">" << fSeq.getName() << " | " << trashCode << endl << contig << endl;
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
        outMisMatch.close();
        outScrapFasta.close();
        if (thisfqualfile != "") {
            inFQual.close();
            inRQual.close();
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
