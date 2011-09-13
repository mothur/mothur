#ifndef CHIMERASLAYERCOMMAND_H
#define CHIMERASLAYERCOMMAND_H

/*
 *  chimeraslayercommand.h
 *  Mothur
 *
 *  Created by westcott on 3/31/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "mothur.h"
#include "command.hpp"
#include "chimera.h"
#include "chimeraslayer.h"

/***********************************************************/

class ChimeraSlayerCommand : public Command {
public:
	ChimeraSlayerCommand(string);
	ChimeraSlayerCommand();
	~ChimeraSlayerCommand() {}
	
	vector<string> setParameters();
	string getCommandName()			{ return "chimera.slayer";		}
	string getCommandCategory()		{ return "Sequence Processing"; }
	string getHelpString();	
	string getCitation() { return "Haas BJ, Gevers D, Earl A, Feldgarden M, Ward DV, Giannokous G, Ciulla D, Tabbaa D, Highlander SK, Sodergren E, Methe B, Desantis TZ, Petrosino JF, Knight R, Birren BW (2011). Chimeric 16S rRNA sequence formation and detection in Sanger and 454-pyrosequenced PCR amplicons. Genome Res. \nhttp://www.mothur.org/wiki/Chimera.slayer"; }
	string getDescription()		{ return "detect chimeric sequences"; }
	
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
	map<string, int> priority;
	
	int driver(linePair*, string, string, string, string);
	int createProcesses(string, string, string, string);
	int divideInHalf(Sequence, string&, string&);
	map<string, int> sortFastaFile(string, string);
		
	#ifdef USE_MPI
	int driverMPI(int, int, MPI_File&, MPI_File&, MPI_File&, MPI_File&, vector<unsigned long long>&);
	#endif

	bool abort, realign, trim, trimera, save;
	string fastafile, templatefile, outputDir, search, namefile, blastlocation;
	int processors, window, iters, increment, numwanted, ksize, match, mismatch, parents, minSimilarity, minCoverage, minBS, minSNP, numSeqs, templateSeqsLength;
	float divR;
	Chimera* chimera;
	
	vector<string> outputNames;
	vector<string> fastaFileNames;
	vector<string> nameFileNames;
	
};

/***********************************************************/

//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
typedef struct slayerData {
	string outputFName; 
	string fasta; 
	string accnos;
	string filename;
	string templatefile;
	string search;
	string blastlocation;
	bool trimera;
	bool trim, realign;
	unsigned long long start;
	unsigned long long end;
	int ksize, match, mismatch, window, minSimilarity, minCoverage, minBS, minSNP, parents, iters, increment, numwanted;
	MothurOut* m;
	float divR;
	map<string, int> priority;
	int count;
	int numNoParents;
	int threadId;
	
	slayerData(){}
	slayerData(string o, string fa, string ac, string f, string te, string se, string bl, bool tri, bool trm, bool re, MothurOut* mout, unsigned long long st, unsigned long long en, int ks, int ma, int mis, int win, int minS, int minC, int miBS, int minSN, int par, int it, int inc, int numw, float div, map<string, int> prior, int tid) {
		outputFName = o;
		fasta = fa;
		accnos = ac;
		filename = f;
		templatefile = te;
		search = se;
		blastlocation = bl;
		trimera = tri;
		trim = trm;
		realign = re;
		m = mout;
		start = st;
		end = en;
		ksize = ks;
		match = ma; 
		mismatch = mis;
		window = win;
		minSimilarity = minS;
		minCoverage = minC;
		minBS = miBS;
		minSNP = minSN;
		parents = par;
		iters = it;
		increment = inc;
		numwanted = numw;
		divR = div;
		priority = prior;
		threadId = tid;
		count = 0;
		numNoParents = 0;
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux)
#else
static DWORD WINAPI MySlayerThreadFunction(LPVOID lpParam){ 
	slayerData* pDataArray;
	pDataArray = (slayerData*)lpParam;
	
	try {
		ofstream out;
		pDataArray->m->openOutputFile(pDataArray->outputFName, out);
		
		ofstream out2;
		pDataArray->m->openOutputFile(pDataArray->accnos, out2);
		
		ofstream out3;
		if (pDataArray->trim) {  pDataArray->m->openOutputFile(pDataArray->fasta, out3); }
		
		ifstream inFASTA;
		pDataArray->m->openInputFile(pDataArray->filename, inFASTA);
		
		
		
		Chimera* chimera;
		if (pDataArray->templatefile != "self") { //you want to run slayer with a reference template
			chimera = new ChimeraSlayer(pDataArray->filename, pDataArray->templatefile, pDataArray->trim, pDataArray->search, pDataArray->ksize, pDataArray->match, pDataArray->mismatch, pDataArray->window, pDataArray->divR, pDataArray->minSimilarity, pDataArray->minCoverage, pDataArray->minBS, pDataArray->minSNP, pDataArray->parents, pDataArray->iters, pDataArray->increment, pDataArray->numwanted, pDataArray->realign, pDataArray->blastlocation, pDataArray->threadId);	
		}else {
			chimera = new ChimeraSlayer(pDataArray->filename, pDataArray->templatefile, pDataArray->trim, pDataArray->priority, pDataArray->search, pDataArray->ksize, pDataArray->match, pDataArray->mismatch, pDataArray->window, pDataArray->divR, pDataArray->minSimilarity, pDataArray->minCoverage, pDataArray->minBS, pDataArray->minSNP, pDataArray->parents, pDataArray->iters, pDataArray->increment, pDataArray->numwanted, pDataArray->realign, pDataArray->blastlocation, pDataArray->threadId);	
		}
		
		//print header if you are process 0
		if ((pDataArray->start == 0) || (pDataArray->start == 1)) {
			chimera->printHeader(out); 
			inFASTA.seekg(0);
		}else { //this accounts for the difference in line endings. 
			inFASTA.seekg(pDataArray->start-1); pDataArray->m->gobble(inFASTA); 
		}
		
		pDataArray->count = pDataArray->end;
		
		if (pDataArray->m->control_pressed) { out.close(); out2.close(); if (pDataArray->trim) { out3.close(); } inFASTA.close(); delete chimera;  return 0;	}
		
		if (chimera->getUnaligned()) { 
			pDataArray->m->mothurOut("Your template sequences are different lengths, please correct."); pDataArray->m->mothurOutEndLine(); 
			out.close(); out2.close(); if (pDataArray->trim) { out3.close(); } inFASTA.close();
			delete chimera;
			return 0; 
		}
		int templateSeqsLength = chimera->getLength();
		
		if (pDataArray->start == 0) { chimera->printHeader(out); }
		
		int count = 0;
		for(int i = 0; i < pDataArray->end; i++){
			
			if (pDataArray->m->control_pressed) {	out.close(); out2.close(); if (pDataArray->trim) { out3.close(); } inFASTA.close(); delete chimera; return 1;	}
			
			Sequence* candidateSeq = new Sequence(inFASTA);  pDataArray->m->gobble(inFASTA);
			string candidateAligned = candidateSeq->getAligned();
			
			if (candidateSeq->getName() != "") { //incase there is a commented sequence at the end of a file
				if (candidateSeq->getAligned().length() != templateSeqsLength) {  
					pDataArray->m->mothurOut(candidateSeq->getName() + " is not the same length as the template sequences. Skipping."); pDataArray->m->mothurOutEndLine();
				}else{
					//find chimeras
					chimera->getChimeras(candidateSeq);
					
					if (pDataArray->m->control_pressed) {	delete candidateSeq; delete chimera; return 1;	}
					
					//if you are not chimeric, then check each half
					data_results wholeResults = chimera->getResults();
					
					//determine if we need to split
					bool isChimeric = false;
					
					if (wholeResults.flag == "yes") {
						string chimeraFlag = "no";
						if(  (wholeResults.results[0].bsa >= pDataArray->minBS && wholeResults.results[0].divr_qla_qrb >= pDataArray->divR)
						   ||
						   (wholeResults.results[0].bsb >= pDataArray->minBS && wholeResults.results[0].divr_qlb_qra >= pDataArray->divR) ) { chimeraFlag = "yes"; }
						
						
						if (chimeraFlag == "yes") {	
							if ((wholeResults.results[0].bsa >= pDataArray->minBS) || (wholeResults.results[0].bsb >= pDataArray->minBS)) { isChimeric = true; }
						}
					}
					
					if ((!isChimeric) && pDataArray->trimera) {
						
						//split sequence in half by bases
						string leftQuery, rightQuery;
						Sequence tempSeq(candidateSeq->getName(), candidateAligned);
						//divideInHalf(tempSeq, leftQuery, rightQuery);
						string queryUnAligned = tempSeq.getUnaligned();
						int numBases = int(queryUnAligned.length() * 0.5);
						
						string queryAligned = tempSeq.getAligned();
						leftQuery = tempSeq.getAligned();
						rightQuery = tempSeq.getAligned();
						
						int baseCount = 0;
						int leftSpot = 0;
						for (int i = 0; i < queryAligned.length(); i++) {
							//if you are a base
							if (isalpha(queryAligned[i])) {		
								baseCount++; 
							}
							
							//if you have half
							if (baseCount >= numBases) {  leftSpot = i; break; } //first half
						}
						
						//blank out right side
						for (int i = leftSpot; i < leftQuery.length(); i++) { leftQuery[i] = '.'; }
						
						//blank out left side
						for (int i = 0; i < leftSpot; i++) { rightQuery[i] = '.'; }
						
						//run chimeraSlayer on each piece
						Sequence* left = new Sequence(candidateSeq->getName(), leftQuery);
						Sequence* right = new Sequence(candidateSeq->getName(), rightQuery);
						
						//find chimeras
						chimera->getChimeras(left);
						data_results leftResults = chimera->getResults();
						
						chimera->getChimeras(right);
						data_results rightResults = chimera->getResults();
						
						//if either piece is chimeric then report
						Sequence trimmed = chimera->print(out, out2, leftResults, rightResults);
						if (pDataArray->trim) { trimmed.printSequence(out3);  }
						
						delete left; delete right;
						
					}else { //already chimeric
						//print results
						Sequence trimmed = chimera->print(out, out2);
						if (pDataArray->trim) { trimmed.printSequence(out3);  }
					}
					
					
				}
				count++;
			}
			
			delete candidateSeq;
			//report progress
			if((count) % 100 == 0){	pDataArray->m->mothurOut("Processing sequence: " + toString(count)); pDataArray->m->mothurOutEndLine();		}
		}
		//report progress
		if((count) % 100 != 0){	pDataArray->m->mothurOut("Processing sequence: " + toString(count)); pDataArray->m->mothurOutEndLine();		}
		
		pDataArray->numNoParents = chimera->getNumNoParents();
		out.close();
		out2.close();
		if (pDataArray->trim) { out3.close(); }
		inFASTA.close();
		delete chimera;
		
		return 0;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "ChimeraSlayerCommand", "MySlayerThreadFunction");
		exit(1);
	}
} 
#endif

/**************************************************************************************************/


#endif


