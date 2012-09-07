#ifndef CLASSIFYSEQSCOMMAND_H
#define CLASSIFYSEQSCOMMAND_H

/*
 *  classifyseqscommand.h
 *  Mothur
 *
 *  Created by westcott on 11/2/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */


#include "command.hpp"
#include "classify.h"
#include "referencedb.h"
#include "sequence.hpp"
#include "bayesian.h"
#include "phylotree.h"
#include "phylosummary.h"
#include "knn.h"


//KNN and Bayesian methods modeled from algorithms in
//Naı¨ve Bayesian Classiﬁer for Rapid Assignment of rRNA Sequences 
//into the New Bacterial Taxonomy􏰎† 
//Qiong Wang,1 George M. Garrity,1,2 James M. Tiedje,1,2 and James R. Cole1* 
//Center for Microbial Ecology1 and Department of Microbiology and Molecular Genetics,2 Michigan State University, 
//East Lansing, Michigan 48824 
//Received 10 January 2007/Accepted 18 June 2007 



class ClassifySeqsCommand : public Command {
	
public:
	ClassifySeqsCommand(string);
	ClassifySeqsCommand();
	~ClassifySeqsCommand();
	
	vector<string> setParameters();
	string getCommandName()			{ return "classify.seqs";		}
	string getCommandCategory()		{ return "Phylotype Analysis";	}
	string getOutputFileNameTag(string, string);
	string getHelpString();	
	string getCitation() { return "Wang Q, Garrity GM, Tiedje JM, Cole JR (2007). Naive Bayesian classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. Appl Environ Microbiol 73: 5261-7. [ for Bayesian classifier ] \nAltschul SF, Madden TL, Schaffer AA, Zhang J, Zhang Z, Miller W, Lipman DJ (1997). Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. Nucleic Acids Res 25: 3389-402. [ for BLAST ] \nDeSantis TZ, Hugenholtz P, Larsen N, Rojas M, Brodie EL, Keller K, Huber T, Dalevi D, Hu P, Andersen GL (2006). Greengenes, a chimera-checked 16S rRNA gene database and workbench compatible with ARB. Appl Environ Microbiol 72: 5069-72. [ for kmer ] \nhttp://www.mothur.org/wiki/Classify.seqs"; }
	string getDescription()		{ return "classify sequences"; }
	
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
	vector<string> fastaFileNames;
	vector<string> namefileNames;
    vector<string> countfileNames;
	vector<string> groupfileNames;
	vector<string> outputNames;
	map<string, vector<string> > nameMap;
	map<string,  vector<string> >::iterator itNames;
	
	Classify* classify;
	ReferenceDB* rdb;
	
	string fastaFileName, templateFileName, countfile, distanceFileName, namefile, search, method, taxonomyFileName, outputDir, groupfile;
	int processors, kmerSize, numWanted, cutoff, iters;
	float match, misMatch, gapOpen, gapExtend;
	bool abort, probs, save, flip, hasName, hasCount;
	
	int driver(linePair*, string, string, string, string);
	int createProcesses(string, string, string, string); 
	string addUnclassifieds(string, int);
	
	int MPIReadNamesFile(string);
	#ifdef USE_MPI
	int driverMPI(int, int, MPI_File&, MPI_File&, MPI_File&, MPI_File&, vector<unsigned long long>&);
	#endif
};

/**************************************************************************************************/
//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct classifyData {
	string taxFName; 
	string tempTFName; 
	string filename;
	string search, taxonomyFileName, templateFileName, method, accnos;
	unsigned long long start;
	unsigned long long end;
	MothurOut* m;
	float match, misMatch, gapOpen, gapExtend;
	int count, kmerSize, threadID, cutoff, iters, numWanted;
	bool probs, flip;
	 
	classifyData(){}
	classifyData(string acc, bool p, string me, string te, string tx, string a, string r, string f, string se, int ks, int i, int numW, MothurOut* mout, unsigned long long st, unsigned long long en, float ma, float misMa, float gapO, float gapE, int cut, int tid, bool fli) {
		accnos = acc;
		taxonomyFileName = tx;
		templateFileName = te;
		taxFName = a;
		tempTFName = r;
		filename = f;
		search = se;
		method = me;
		m = mout;
		start = st;
		end = en;
		match = ma; 
		misMatch = misMa;
		gapOpen = gapO; 
		gapExtend = gapE; 
		kmerSize = ks;
		cutoff = cut;
		iters = i;
		numWanted = numW;
		threadID = tid;
		probs = p;
		count = 0;
		flip = fli;
	}
};

/**************************************************************************************************/
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
static DWORD WINAPI MyClassThreadFunction(LPVOID lpParam){ 
	classifyData* pDataArray;
	pDataArray = (classifyData*)lpParam;
	
	try {
		ofstream outTax;
		pDataArray->m->openOutputFile(pDataArray->taxFName, outTax);
		
		ofstream outTaxSimple;
		pDataArray->m->openOutputFile(pDataArray->tempTFName, outTaxSimple);
		
		ofstream outAcc;
		pDataArray->m->openOutputFile(pDataArray->accnos, outAcc);
		
		ifstream inFASTA;
		pDataArray->m->openInputFile(pDataArray->filename, inFASTA);
		
		string taxonomy;
				
		//print header if you are process 0
		if ((pDataArray->start == 0) || (pDataArray->start == 1)) {
			inFASTA.seekg(0);
		}else { //this accounts for the difference in line endings. 
			inFASTA.seekg(pDataArray->start-1); pDataArray->m->gobble(inFASTA); 
		}
		
		pDataArray->count = pDataArray->end;
		
		//make classify
		Classify* myclassify;
		if(pDataArray->method == "bayesian"){	myclassify = new Bayesian(pDataArray->taxonomyFileName, pDataArray->templateFileName, pDataArray->search, pDataArray->kmerSize, pDataArray->cutoff, pDataArray->iters, pDataArray->threadID, pDataArray->flip);		}
		else if(pDataArray->method == "knn"){	myclassify = new Knn(pDataArray->taxonomyFileName, pDataArray->templateFileName, pDataArray->search, pDataArray->kmerSize, pDataArray->gapOpen, pDataArray->gapExtend, pDataArray->match, pDataArray->misMatch, pDataArray->numWanted, pDataArray->threadID);				}
		else {
			pDataArray->m->mothurOut(pDataArray->search + " is not a valid method option. I will run the command using bayesian.");
			pDataArray->m->mothurOutEndLine();
			myclassify = new Bayesian(pDataArray->taxonomyFileName, pDataArray->templateFileName, pDataArray->search, pDataArray->kmerSize, pDataArray->cutoff, pDataArray->iters, pDataArray->threadID, pDataArray->flip);	
		}
		
		if (pDataArray->m->control_pressed) { delete myclassify; return 0; }
		
		int count = 0;
		for(int i = 0; i < pDataArray->end; i++){ //end is the number of sequences to process
			
			if (pDataArray->m->control_pressed) { delete myclassify; return 0; }
			
			Sequence* candidateSeq = new Sequence(inFASTA); pDataArray->m->gobble(inFASTA);
			
			if (candidateSeq->getName() != "") {
				
				taxonomy = myclassify->getTaxonomy(candidateSeq);
				
				if (pDataArray->m->control_pressed) { delete candidateSeq; return 0; }
				
				if (taxonomy == "unknown;") { pDataArray->m->mothurOut("[WARNING]: " + candidateSeq->getName() + " could not be classified. You can use the remove.lineage command with taxon=unknown; to remove such sequences."); pDataArray->m->mothurOutEndLine(); }

				//output confidence scores or not
				if (pDataArray->probs) {
					outTax << candidateSeq->getName() << '\t' << taxonomy << endl;
				}else{
					outTax << candidateSeq->getName() << '\t' << myclassify->getSimpleTax() << endl;
				}
					
				outTaxSimple << candidateSeq->getName() << '\t' << myclassify->getSimpleTax() << endl;
					
				if (myclassify->getFlipped()) { outAcc << candidateSeq->getName() << endl; }
				
				count++;
			}
			delete candidateSeq;
			//report progress
			if((count) % 100 == 0){	pDataArray->m->mothurOut("Processing sequence: " + toString(count)); pDataArray->m->mothurOutEndLine();		}
			
		}
		//report progress
		if((count) % 100 != 0){	pDataArray->m->mothurOut("Processing sequence: " + toString(count)); pDataArray->m->mothurOutEndLine();		}
		
		delete myclassify;
		inFASTA.close();
		outTax.close();
		outTaxSimple.close();
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "ClassifySeqsCommand", "MyClassThreadFunction");
		exit(1);
	}
} 
#endif




#endif

