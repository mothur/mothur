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


struct fastqRead {
	vector<int> scores;
	string name;
	string sequence;
	
	fastqRead() { name = ""; sequence = ""; scores.clear(); };
	fastqRead(string n, string s, vector<int> sc) : name(n), sequence(s), scores(sc) {};
	~fastqRead() {};
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
    string getOutputFileNameTag(string, string);
	string getHelpString();	
    string getCitation() { return "http://www.mothur.org/wiki/Make.contigs"; }
    string getDescription()		{ return "description"; }
    
    int execute(); 
    void help() { m->mothurOut(getHelpString()); }	
    
private:
    bool abort;
    string outputDir, ffastqfile, rfastqfile, align;
	float match, misMatch, gapOpen, gapExtend;
	int processors, longestBase, threshold;
    vector<string> outputNames;
    
    fastqRead readFastq(ifstream&);
    vector< vector<string> > readFastqFiles(int&);
    bool checkReads(fastqRead&, fastqRead&);
    int createProcesses(vector< vector<string> >, string, string, string);
    int driver(vector<string>, string, string, string);
};

/**************************************************************************************************/

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct contigsData {
	string outputFasta; 
	string outputQual; 
	string outputMisMatches;
	string align;
    vector<string> files;
	MothurOut* m;
	float match, misMatch, gapOpen, gapExtend;
	int count, threshold, threadID;
	
	contigsData(){}
	contigsData(vector<string> f, string of, string oq, string om, string al, MothurOut* mout, float ma, float misMa, float gapO, float gapE, int thr, int tid) {
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
        
        ifstream inFFasta, inRFasta, inFQual, inRQual;
        pDataArray->m->openInputFile(thisffastafile, inFFasta);
        pDataArray->m->openInputFile(thisfqualfile, inFQual);
        pDataArray->m->openInputFile(thisrfastafile, inRFasta);
        pDataArray->m->openInputFile(thisrqualfile, inRQual);
        
        ofstream outFasta, outQual, outMisMatch;
        pDataArray->m->openOutputFile(pDataArray->outputFasta, outFasta);
        pDataArray->m->openOutputFile(pDataArray->outputQual, outQual);
        pDataArray->m->openOutputFile(pDataArray->outputMisMatches, outMisMatch);
        outMisMatch << "Name\tLength\tMisMatches\n";
        
        while ((!inFQual.eof()) && (!inFFasta.eof()) && (!inRFasta.eof()) && (!inRQual.eof())) {
            
            if (pDataArray->m->control_pressed) { break; }
            
            //read seqs and quality info
            Sequence fSeq(inFFasta); pDataArray->m->gobble(inFFasta);
            Sequence rSeq(inRFasta); pDataArray->m->gobble(inRFasta);
            QualityScores fQual(inFQual); pDataArray->m->gobble(inFQual);
            QualityScores rQual(inRQual); pDataArray->m->gobble(inRQual);
            
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
            
            for (int i = 0; i < length; i++) {
                if (seq1[i] == seq2[i]) { //match, add base and choose highest score
                    contig += seq1[i];
                    contigScores.push_back(scores1[ABaseMap[i]]);
                    if (scores1[ABaseMap[i]] < scores2[BBaseMap[i]]) { contigScores[i] = scores2[BBaseMap[i]]; }
                }else if (((seq1[i] == '.') || (seq1[i] == '-')) && ((seq2[i] != '-') && (seq2[i] != '.'))) { //seq1 is a gap and seq2 is a base, choose seq2, unless quality score for base is below threshold. In that case eliminate base
                    if (scores2[BBaseMap[i]] >= pDataArray->threshold)) {
                        contig += seq2[i];
                        contigScores.push_back(scores2[BBaseMap[i]]);
                    }
                }else if (((seq2[i] == '.') || (seq2[i] == '-')) && ((seq1[i] != '-') && (seq1[i] != '.'))) { //seq2 is a gap and seq1 is a base, choose seq1, unless quality score for base is below threshold. In that case eliminate base
                    if (scores1[ABaseMap[i]] >= pDataArray->threshold) { 
                        contig += seq1[i];
                        contigScores.push_back(scores1[ABaseMap[i]]);
                    }
                }else if (((seq1[i] != '-') && (seq1[i] != '.')) && ((seq2[i] != '-') && (seq2[i] != '.'))) { //both bases choose one with better quality
                    char c = seq1[i];
                    contigScores.push_back(scores1[ABaseMap[i]]);
                    if (scores1[ABaseMap[i]] < scores2[BBaseMap[i]]) { contigScores[i] = scores2[BBaseMap[i]]; c = seq2[i]; }
                    contig += c;
                    numMismatches++;
                }else { //should never get here
                    pDataArray->m->mothurOut("[ERROR]: case I didn't think of seq1 = " + toString(seq1[i]) + " and seq2 = " + toString(seq2[i]) + "\n");
                }
            }
            
            //output
            outFasta << ">" << fSeq.getName() << endl << contig << endl;
            outQual << ">" << fSeq.getName() << endl;
            for (int i = 0; i < contigScores.size(); i++) { outQual << contigScores[i] << ' '; }
            outQual << endl;
            outMisMatch << fSeq.getName() << '\t' << contig.length() << '\t' << numMismatches << endl;
            
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
        delete alignment;
        
        if (pDataArray->m->control_pressed) { pDataArray->m->mothurRemove(pDataArray->outputQual); pDataArray->m->mothurRemove(pDataArray->outputFasta);  pDataArray->m->mothurRemove(pDataArray->outputMisMatches);}
        
        return 0;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "AlignCommand", "MyContigsThreadFunction");
		exit(1);
	}
} 
#endif


#endif
