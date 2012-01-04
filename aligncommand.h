#ifndef ALIGNCOMMAND_H
#define ALIGNCOMMAND_H

/*
 *  aligncommand.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 5/15/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "database.hpp"
#include "alignment.hpp"
#include "alignmentdb.h"
#include "sequence.hpp"

#include "gotohoverlap.hpp"
#include "needlemanoverlap.hpp"
#include "blastalign.hpp"
#include "noalign.hpp"

#include "nast.hpp"
#include "nastreport.hpp"


class AlignCommand : public Command {
	
public:
	AlignCommand(string);	
	AlignCommand();
	~AlignCommand();
	
	vector<string> setParameters();
	string getCommandName()			{ return "align.seqs";			}
	string getCommandCategory()		{ return "Sequence Processing"; }
	string getHelpString();	
	string getCitation() { return "DeSantis TZ, Jr., Hugenholtz P, Keller K, Brodie EL, Larsen N, Piceno YM, Phan R, Andersen GL (2006). NAST: a multiple sequence alignment server for comparative analysis of 16S rRNA genes. Nucleic Acids Res 34: W394-9.\nSchloss PD (2009). A high-throughput DNA sequence aligner for microbial ecology studies. PLoS ONE 4: e8230.\nSchloss PD (2010). The effects of alignment quality, distance calculation method, sequence filtering, and region on the analysis of 16S rRNA gene-based studies. PLoS Comput Biol 6: e1000844.\nhttp://www.mothur.org/wiki/Align.seqs http://www.mothur.org/wiki/Align.seqs"; }
	string getDescription()		{ return "align sequences"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }	
	
private:
	struct linePair {
		unsigned long long start;
		unsigned long long end;
		linePair(unsigned long long i, unsigned long long j) : start(i), end(j) {}
	};
	vector<int> processIDS;   //processid
	vector<linePair*> lines;
	bool MPIWroteAccnos;
	
	AlignmentDB* templateDB;
	
	int driver(linePair*, string, string, string, string);
	int createProcesses(string, string, string, string);
	void appendAlignFiles(string, string); 
	void appendReportFiles(string, string);
	
	#ifdef USE_MPI
	int driverMPI(int, int, MPI_File&, MPI_File&, MPI_File&, MPI_File&, vector<unsigned long long>&);
	#endif
	
	string candidateFileName, templateFileName, distanceFileName, search, align, outputDir;
	float match, misMatch, gapOpen, gapExtend, threshold;
	int processors, kmerSize;
	vector<string> candidateFileNames;
	vector<string> outputNames;
	
	bool abort, flip, calledHelp, save;

};

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct alignData {
	string templateFileName;
	string alignFName; 
	string reportFName; 
	string accnosFName;
	string filename;
	string align;
	string search;
	bool flip;
	unsigned long long start;
	unsigned long long end;
	MothurOut* m;
	//AlignmentDB* templateDB;
	float match, misMatch, gapOpen, gapExtend, threshold;
	int count, kmerSize, threadID;
	
	alignData(){}
	alignData(string te, string a, string r, string ac, string f, string al, string se, int ks, MothurOut* mout, unsigned long long st, unsigned long long en, bool fl, float ma, float misMa, float gapO, float gapE, float thr, int tid) {
		templateFileName = te;
		alignFName = a;
		reportFName = r;
		accnosFName = ac;
		filename = f;
		flip = fl;
		m = mout;
		start = st;
		end = en;
		//templateDB = tdb;
		match = ma; 
		misMatch = misMa;
		gapOpen = gapO; 
		gapExtend = gapE; 
		threshold = thr;
		align = al;
		search = se;
		count = 0;
		kmerSize = ks;
		threadID = tid;
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
#else
static DWORD WINAPI MyAlignThreadFunction(LPVOID lpParam){ 
	alignData* pDataArray;
	pDataArray = (alignData*)lpParam;
	
	try {
		ofstream alignmentFile;
		pDataArray->m->openOutputFile(pDataArray->alignFName, alignmentFile);
		
		ofstream accnosFile;
		pDataArray->m->openOutputFile(pDataArray->accnosFName, accnosFile);
		
		NastReport report(pDataArray->reportFName);
		
		ifstream inFASTA;
		pDataArray->m->openInputFile(pDataArray->filename, inFASTA);
		
		//print header if you are process 0
		if ((pDataArray->start == 0) || (pDataArray->start == 1)) {
			inFASTA.seekg(0);
		}else { //this accounts for the difference in line endings. 
			inFASTA.seekg(pDataArray->start-1); pDataArray->m->gobble(inFASTA); 
		}
		
		pDataArray->count = pDataArray->end;
		
		AlignmentDB* templateDB = new AlignmentDB(templateFileName, pDataArray->search, pDataArray->kmerSize, pDataArray->gapOpen, pDataArray->gapExtend, pDataArray->match, pDataArray->misMatch, pDataArray->threadID);
		
		//moved this into driver to avoid deep copies in windows paralellized version
		Alignment* alignment;
		int longestBase = templateDB->getLongestBase();
		if(pDataArray->align == "gotoh")			{	alignment = new GotohOverlap(pDataArray->gapOpen, pDataArray->gapExtend, pDataArray->match, pDataArray->misMatch, longestBase);			}
		else if(pDataArray->align == "needleman")	{	alignment = new NeedlemanOverlap(pDataArray->gapOpen, pDataArray->match, pDataArray->misMatch, longestBase);				}
		else if(pDataArray->align == "blast")		{	alignment = new BlastAlignment(pDataArray->gapOpen, pDataArray->gapExtend, pDataArray->match, pDataArray->misMatch);		}
		else if(pDataArray->align == "noalign")		{	alignment = new NoAlign();													}
		else {
			pDataArray->m->mothurOut(pDataArray->align + " is not a valid alignment option. I will run the command using needleman.");
			pDataArray->m->mothurOutEndLine();
			alignment = new NeedlemanOverlap(pDataArray->gapOpen, pDataArray->match, pDataArray->misMatch, longestBase);
		}
		
		int count = 0;
		for(int i = 0; i < pDataArray->end; i++){ //end is the number of sequences to process
			
			if (pDataArray->m->control_pressed) {  break; }
			
			Sequence* candidateSeq = new Sequence(inFASTA);  pDataArray->m->gobble(inFASTA);
			report.setCandidate(candidateSeq);
			
			int origNumBases = candidateSeq->getNumBases();
			string originalUnaligned = candidateSeq->getUnaligned();
			int numBasesNeeded = origNumBases * pDataArray->threshold;
			
			if (candidateSeq->getName() != "") { //incase there is a commented sequence at the end of a file
				if (candidateSeq->getUnaligned().length() > alignment->getnRows()) {
					alignment->resize(candidateSeq->getUnaligned().length()+1);
				}
				
				Sequence temp = templateDB->findClosestSequence(candidateSeq);
				Sequence* templateSeq = &temp;
				
				float searchScore = templateDB->getSearchScore();
				
				Nast* nast = new Nast(alignment, candidateSeq, templateSeq);
				
				Sequence* copy;
				
				Nast* nast2;
				bool needToDeleteCopy = false;  //this is needed in case you have you enter the ifs below
				//since nast does not make a copy of hte sequence passed, and it is used by the reporter below
				//you can't delete the copy sequence til after you report, but you may choose not to create it in the first place
				//so this bool tells you if you need to delete it
				
				//if there is a possibility that this sequence should be reversed
				if (candidateSeq->getNumBases() < numBasesNeeded) {
					
					string wasBetter =  "";
					//if the user wants you to try the reverse
					if (pDataArray->flip) {
						
						//get reverse compliment
						copy = new Sequence(candidateSeq->getName(), originalUnaligned);
						copy->reverseComplement();
						
						//rerun alignment
						Sequence temp2 = templateDB->findClosestSequence(copy);
						Sequence* templateSeq2 = &temp2;
						
						searchScore = templateDB->getSearchScore();
						
						nast2 = new Nast(alignment, copy, templateSeq2);
						
						//check if any better
						if (copy->getNumBases() > candidateSeq->getNumBases()) {
							candidateSeq->setAligned(copy->getAligned());  //use reverse compliments alignment since its better
							templateSeq = templateSeq2; 
							delete nast;
							nast = nast2;
							needToDeleteCopy = true;
							wasBetter = "\treverse complement produced a better alignment, so mothur used the reverse complement.";
						}else{  
							wasBetter = "\treverse complement did NOT produce a better alignment so it was not used, please check sequence.";
							delete nast2;
							delete copy;	
						}
					}
					
					//create accnos file with names
					accnosFile << candidateSeq->getName() << wasBetter << endl;
				}
				
				report.setTemplate(templateSeq);
				report.setSearchParameters(pDataArray->search, searchScore);
				report.setAlignmentParameters(pDataArray->align, alignment);
				report.setNastParameters(*nast);
				
				alignmentFile << '>' << candidateSeq->getName() << '\n' << candidateSeq->getAligned() << endl;
				
				report.print();
				delete nast;
				if (needToDeleteCopy) {   delete copy;   }
				
				count++;
			}
			delete candidateSeq;
			
			//report progress
			if((count) % 100 == 0){	pDataArray->m->mothurOut(toString(count)); pDataArray->m->mothurOutEndLine();		}
			
		}
		//report progress
		if((count) % 100 != 0){	pDataArray->m->mothurOut(toString(count)); pDataArray->m->mothurOutEndLine();		}
		
		delete alignment;
		delete templateDB;
		alignmentFile.close();
		inFASTA.close();
		accnosFile.close();
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "AlignCommand", "MyAlignThreadFunction");
		exit(1);
	}
} 
#endif



#endif
