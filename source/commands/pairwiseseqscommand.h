#ifndef PAIRWISESEQSCOMMAND_H
#define PAIRWISESEQSCOMMAND_H

/*
 *  pairwiseseqscommand.h
 *  Mothur
 *
 *  Created by westcott on 10/20/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "database.hpp"
#include "alignment.hpp"
#include "validcalculator.h"
#include "dist.h"
#include "sequencedb.h"
#include "sequence.hpp"

#include "gotohoverlap.hpp"
#include "needlemanoverlap.hpp"
#include "blastalign.hpp"
#include "noalign.hpp"

#include "ignoregaps.h"
#include "eachgapdist.h"
#include "eachgapignore.h"
#include "onegapdist.h"
#include "onegapignore.h"

class PairwiseSeqsCommand : public Command {
	
public:
	PairwiseSeqsCommand(string);	
	PairwiseSeqsCommand();
	~PairwiseSeqsCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "pairwise.seqs";		}
	string getCommandCategory()		{ return "Sequence Processing"; }
	
	string getHelpString();	
    string getOutputPattern(string);	
	string getCitation() { return "Needleman SB, Wunsch CD (1970). A general method applicable to the search for similarities in the amino acid sequence of two proteins. J Mol Biol 48: 443-53. [ for needleman ]\nGotoh O (1982). An improved algorithm for matching biological sequences. J Mol Biol 162: 705-8. [ for gotoh ] \nhttp://www.mothur.org/wiki/Pairwise.seqs"; }
	string getDescription()		{ return "calculates pairwise distances from an unaligned fasta file"; }

	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	struct distlinePair {
		int start;
		int end;
	};
	
	vector<int> processIDS;   //end line, processid
	vector<distlinePair> lines;
	
	SequenceDB alignDB;
	
	void createProcesses(string);
	int driver(int, int, string, float);
	int driver(int, int, string, string);
	
	string fastaFileName, align, calc, outputDir, output;
	float match, misMatch, gapOpen, gapExtend, cutoff;
	int processors, longestBase;
	vector<string> fastaFileNames, Estimators;
	vector<string> outputNames;
	
	bool abort, countends, compress;
};

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct pairwiseData {
    string outputFileName;
	string align, square, distcalcType, output;
	unsigned long long start;
	unsigned long long end;
	MothurOut* m;
	float match, misMatch, gapOpen, gapExtend, cutoff;
	int count, threadID, longestBase;
    bool countends;
    SequenceDB alignDB;
	
	pairwiseData(){}
	pairwiseData(string ofn, string al, string sq, string di, bool co, string op, SequenceDB DB, MothurOut* mout, unsigned long long st, unsigned long long en, float ma, float misMa, float gapO, float gapE, int thr, float cu, int tid) {
		outputFileName = ofn;
		m = mout;
		start = st;
		end = en;
		match = ma; 
		misMatch = misMa;
		gapOpen = gapO; 
		gapExtend = gapE; 
		longestBase = thr;
		align = al;
        square = sq;
        distcalcType = di;
        countends = co;
        alignDB = DB;
		count = 0;
        output = op;
        cutoff = cu;
		threadID = tid;
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyPairwiseSquareThreadFunction(LPVOID lpParam){ 
	pairwiseData* pDataArray;
	pDataArray = (pairwiseData*)lpParam;
	
	try {
		ofstream outFile((pDataArray->outputFileName).c_str(), ios::trunc);
		outFile.setf(ios::fixed, ios::showpoint);
		outFile << setprecision(4);
		
		pDataArray->count = 0;
        
        int startTime = time(NULL);
        
        Alignment* alignment;
        if(pDataArray->align == "gotoh")			{	alignment = new GotohOverlap(pDataArray->gapOpen, pDataArray->gapExtend, pDataArray->match, pDataArray->misMatch, pDataArray->longestBase);			}
		else if(pDataArray->align == "needleman")	{	alignment = new NeedlemanOverlap(pDataArray->gapOpen, pDataArray->match, pDataArray->misMatch, pDataArray->longestBase);				}
		else if(pDataArray->align == "blast")		{	alignment = new BlastAlignment(pDataArray->gapOpen, pDataArray->gapExtend, pDataArray->match, pDataArray->misMatch);		}
		else if(pDataArray->align == "noalign")		{	alignment = new NoAlign();													}
		else {
			pDataArray->m->mothurOut(pDataArray->align + " is not a valid alignment option. I will run the command using needleman.");
			pDataArray->m->mothurOutEndLine();
			alignment = new NeedlemanOverlap(pDataArray->gapOpen, pDataArray->match, pDataArray->misMatch, pDataArray->longestBase);
		}
		
        ValidCalculators validCalculator;
        Dist* distCalculator;
        if (pDataArray->countends) {
            if (validCalculator.isValidCalculator("distance", pDataArray->distcalcType) == true) { 
                if (pDataArray->distcalcType == "nogaps")			{	distCalculator = new ignoreGaps();	}
                else if (pDataArray->distcalcType == "eachgap")	{	distCalculator = new eachGapDist();	}
                else if (pDataArray->distcalcType == "onegap")		{	distCalculator = new oneGapDist();	}
            }
        }else {
            if (validCalculator.isValidCalculator("distance", pDataArray->distcalcType) == true) { 
                if (pDataArray->distcalcType == "nogaps")		{	distCalculator = new ignoreGaps();					}
                else if (pDataArray->distcalcType == "eachgap"){	distCalculator = new eachGapIgnoreTermGapDist();	}
                else if (pDataArray->distcalcType == "onegap")	{	distCalculator = new oneGapIgnoreTermGapDist();		}
            }
        }

        if(pDataArray->start == 0){	outFile << pDataArray->alignDB.getNumSeqs() << endl;	}
		
		for(int i=pDataArray->start;i<pDataArray->end;i++){
            pDataArray->count++;
            
			string name = pDataArray->alignDB.get(i).getName();
			//pad with spaces to make compatible
			if (name.length() < 10) { while (name.length() < 10) {  name += " ";  } }
            
			outFile << name;
			
			for(int j=0;j<pDataArray->alignDB.getNumSeqs();j++){
				
				if (pDataArray->m->control_pressed) { outFile.close(); delete alignment; delete distCalculator; return 0;  }
				
				if (pDataArray->alignDB.get(i).getUnaligned().length() > alignment->getnRows()) {
					alignment->resize(pDataArray->alignDB.get(i).getUnaligned().length()+1);
				}
				
				if (pDataArray->alignDB.get(j).getUnaligned().length() > alignment->getnRows()) {
					alignment->resize(pDataArray->alignDB.get(j).getUnaligned().length()+1);
				}
				
				Sequence seqI(pDataArray->alignDB.get(i).getName(), pDataArray->alignDB.get(i).getAligned());
				Sequence seqJ(pDataArray->alignDB.get(j).getName(), pDataArray->alignDB.get(j).getAligned());
				
				alignment->align(seqI.getUnaligned(), seqJ.getUnaligned());
				seqI.setAligned(alignment->getSeqAAln());
				seqJ.setAligned(alignment->getSeqBAln());
				
				distCalculator->calcDist(seqI, seqJ);
				double dist = distCalculator->getDist();
                
                if (pDataArray->m->debug) { pDataArray->m->mothurOut("[DEBUG]: " + seqI.getName() + '\t' +  alignment->getSeqAAln() + '\n' + seqJ.getName() + alignment->getSeqBAln() + '\n' + "distance = " + toString(dist) + "\n"); }
                
				outFile  << '\t' << dist;
			}
			
			outFile << endl; 
			
			if(i % 100 == 0){
				pDataArray->m->mothurOutJustToScreen(toString(i) + "\t" + toString(time(NULL) - startTime)+"\n"); 
			}
			
		}
		pDataArray->m->mothurOutJustToScreen(toString(pDataArray->count) + "\t" + toString(time(NULL) - startTime)+"\n");
		
		outFile.close();
        delete alignment;
        delete distCalculator;

        
    }
	catch(exception& e) {
		pDataArray->m->errorOut(e, "PairwiseSeqsCommand", "MyPairwiseSquareThreadFunction");
		exit(1);
	}
} 

/**************************************************************************************************/
static DWORD WINAPI MyPairwiseThreadFunction(LPVOID lpParam){ 
	pairwiseData* pDataArray;
	pDataArray = (pairwiseData*)lpParam;
	
	try {
		ofstream outFile((pDataArray->outputFileName).c_str(), ios::trunc);
		outFile.setf(ios::fixed, ios::showpoint);
		outFile << setprecision(4);
        
        int startTime = time(NULL);
        
        Alignment* alignment;
        if(pDataArray->align == "gotoh")			{	alignment = new GotohOverlap(pDataArray->gapOpen, pDataArray->gapExtend, pDataArray->match, pDataArray->misMatch, pDataArray->longestBase);			}
		else if(pDataArray->align == "needleman")	{	alignment = new NeedlemanOverlap(pDataArray->gapOpen, pDataArray->match, pDataArray->misMatch, pDataArray->longestBase);				}
		else if(pDataArray->align == "blast")		{	alignment = new BlastAlignment(pDataArray->gapOpen, pDataArray->gapExtend, pDataArray->match, pDataArray->misMatch);		}
		else if(pDataArray->align == "noalign")		{	alignment = new NoAlign();													}
		else {
			pDataArray->m->mothurOut(pDataArray->align + " is not a valid alignment option. I will run the command using needleman.");
			pDataArray->m->mothurOutEndLine();
			alignment = new NeedlemanOverlap(pDataArray->gapOpen, pDataArray->match, pDataArray->misMatch, pDataArray->longestBase);
		}
		
        ValidCalculators validCalculator;
        Dist* distCalculator;
        if (pDataArray->countends) {
            if (validCalculator.isValidCalculator("distance", pDataArray->distcalcType) == true) { 
                if (pDataArray->distcalcType == "nogaps")			{	distCalculator = new ignoreGaps();	}
                else if (pDataArray->distcalcType == "eachgap")	{	distCalculator = new eachGapDist();	}
                else if (pDataArray->distcalcType == "onegap")		{	distCalculator = new oneGapDist();	}
            }
        }else {
            if (validCalculator.isValidCalculator("distance", pDataArray->distcalcType) == true) { 
                if (pDataArray->distcalcType == "nogaps")		{	distCalculator = new ignoreGaps();					}
                else if (pDataArray->distcalcType == "eachgap"){	distCalculator = new eachGapIgnoreTermGapDist();	}
                else if (pDataArray->distcalcType == "onegap")	{	distCalculator = new oneGapIgnoreTermGapDist();		}
            }
        }
        
        if((pDataArray->output == "lt") && pDataArray->start == 0){	outFile << pDataArray->alignDB.getNumSeqs() << endl;	}
		
        pDataArray->count = 0;
		for(int i=pDataArray->start;i<pDataArray->end;i++){
            pDataArray->count++;
            
			if(pDataArray->output == "lt")	{	
				string name = pDataArray->alignDB.get(i).getName();
				if (name.length() < 10) { //pad with spaces to make compatible
					while (name.length() < 10) {  name += " ";  }
				}
				outFile << name;
			}

			
			for(int j=0;j<i;j++){
				
				if (pDataArray->m->control_pressed) { outFile.close(); delete alignment; delete distCalculator; return 0;  }
				
				if (pDataArray->alignDB.get(i).getUnaligned().length() > alignment->getnRows()) {
					alignment->resize(pDataArray->alignDB.get(i).getUnaligned().length()+1);
				}
				
				if (pDataArray->alignDB.get(j).getUnaligned().length() > alignment->getnRows()) {
					alignment->resize(pDataArray->alignDB.get(j).getUnaligned().length()+1);
				}
				
				Sequence seqI(pDataArray->alignDB.get(i).getName(), pDataArray->alignDB.get(i).getAligned());
				Sequence seqJ(pDataArray->alignDB.get(j).getName(), pDataArray->alignDB.get(j).getAligned());
				
				alignment->align(seqI.getUnaligned(), seqJ.getUnaligned());
				seqI.setAligned(alignment->getSeqAAln());
				seqJ.setAligned(alignment->getSeqBAln());
				
				distCalculator->calcDist(seqI, seqJ);
				double dist = distCalculator->getDist();
                
                if (pDataArray->m->debug) { pDataArray->m->mothurOut("[DEBUG]: " + seqI.getName() + '\t' +  alignment->getSeqAAln() + '\n' + seqJ.getName() + alignment->getSeqBAln() + '\n' + "distance = " + toString(dist) + "\n"); }
                
				if(dist <= pDataArray->cutoff){
					if (pDataArray->output == "column") { outFile << pDataArray->alignDB.get(i).getName() << ' ' << pDataArray->alignDB.get(j).getName() << ' ' << dist << endl; }
				}
				if (pDataArray->output == "lt") {  outFile << '\t' << dist; }
			}
			
			if (pDataArray->output == "lt") { outFile << endl; }
			
			if(i % 100 == 0){
				pDataArray->m->mothurOutJustToScreen(toString(i) + "\t" + toString(time(NULL) - startTime)+"\n"); 
			}
			
		}
		pDataArray->m->mothurOutJustToScreen(toString(pDataArray->end-1) + "\t" + toString(time(NULL) - startTime)+"\n");
		
		outFile.close();
        delete alignment;
        delete distCalculator;
        
        
    }
	catch(exception& e) {
		pDataArray->m->errorOut(e, "PairwiseSeqsCommand", "MyPairwiseThreadFunction");
		exit(1);
	}
} 

#endif


#endif

