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
#include "sequenceparser.h"
#include "sequencecountparser.h"

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
    string getOutputPattern(string);	
	string getCitation() { return "Haas BJ, Gevers D, Earl A, Feldgarden M, Ward DV, Giannokous G, Ciulla D, Tabbaa D, Highlander SK, Sodergren E, Methe B, Desantis TZ, Petrosino JF, Knight R, Birren BW (2011). Chimeric 16S rRNA sequence formation and detection in Sanger and 454-pyrosequenced PCR amplicons. Genome Res  21:494.\nhttp://www.mothur.org/wiki/Chimera.slayer"; }
	string getDescription()		{ return "detect chimeric sequences"; }
	
	int execute(); 
	void help() { m->mothurOut(getHelpString()); }		
	
private:

	vector<int> processIDS;   //processid
	vector<linePair> lines;
	
	int driver(linePair, string, string, string, string, map<string, int>&);
	int createProcesses(string, string, string, string, map<string, int>&);
	int divideInHalf(Sequence, string&, string&);
	map<string, int> sortFastaFile(string, string);
	map<string, int> sortFastaFile(vector<Sequence>&, map<string, string>&, string newFile);
    int sortFastaFile(vector<Sequence>&, map<string, int>&, string newFile);
	string getNamesFile(string&);
	//int setupChimera(string,);
	int MPIExecute(string, string, string, string, map<string, int>&);
	int deconvoluteResults(map<string, string>&, string, string, string);
	map<string, int> priority;
	int setUpForSelfReference(SequenceParser*&, map<string, string>&, map<string, map<string, int> >&, int);
    int setUpForSelfReference(SequenceCountParser*&, map<string, string>&, map<string, map<string, int> >&, int);
	int driverGroups(string, string, string, map<string, map<string, int> >&, map<string, string>&, string);
	int createProcessesGroups(string, string, string, map<string, map<string, int> >&, map<string, string>&, string, string);
	int MPIExecuteGroups(string, string, string, map<string, map<string, int> >&, map<string, string>&, string, string);

		
	#ifdef USE_MPI
	int driverMPI(int, int, MPI_File&, MPI_File&, MPI_File&, MPI_File&, set<string>&, vector<unsigned long long>&, string, map<string, int>&, bool);
	#endif

	bool abort, realign, trim, trimera, save, hasName, hasCount, dups;
	string fastafile, groupfile, templatefile, outputDir, search, namefile, countfile, blastlocation;
	int processors, window, iters, increment, numwanted, ksize, match, mismatch, parents, minSimilarity, minCoverage, minBS, minSNP, numSeqs, templateSeqsLength;
	float divR;
	
    map<string, map<string, string> > group2NameMap;
	vector<string> outputNames;
	vector<string> fastaFileNames;
	vector<string> nameFileNames;
	vector<string> groupFileNames;
	
};

/***********************************************************/

//custom data structure for threads to use.
// This is passed by void pointer so it can be any data type
// that can be passed using a single void pointer (LPVOID).
struct slayerData {
	string outputFName; 
	string fasta; 
	string accnos;
	string filename, countlist;
	string templatefile;
	string search;
	string blastlocation;
	bool trimera;
	bool trim, realign, dups, hasCount;
	unsigned long long start;
	unsigned long long end;
	int ksize, match, mismatch, window, minSimilarity, minCoverage, minBS, minSNP, parents, iters, increment, numwanted;
	MothurOut* m;
	float divR;
	map<string, int> priority;
	int count;
	int numNoParents;
	int threadId;
	map<string, map<string, int> > fileToPriority;
	map<string, string> fileGroup;
    map<string, map<string, string> > group2NameMap;
	
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
	slayerData(map<string, map<string, string> > g2n, bool hc, bool dps, string cl, string o, string fa, string ac, string te, string se, string bl, bool tri, bool trm, bool re, MothurOut* mout, map<string, map<string, int> >& fPriority, map<string, string>& fileG, int ks, int ma, int mis, int win, int minS, int minC, int miBS, int minSN, int par, int it, int inc, int numw, float div, map<string, int> prior, int tid) {
		outputFName = o;
		fasta = fa;
		accnos = ac;
		templatefile = te;
		search = se;
		blastlocation = bl;
        countlist = cl;
        dups = dps;
        hasCount = hc;
        group2NameMap = g2n;
		trimera = tri;
		trim = trm;
		realign = re;
		m = mout;
		fileGroup = fileG;
		fileToPriority = fPriority;
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
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
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
            pDataArray->m->zapGremlins(inFASTA);
		}else { //this accounts for the difference in line endings. 
			inFASTA.seekg(pDataArray->start-1); pDataArray->m->gobble(inFASTA); 
		}
		
		if (pDataArray->m->control_pressed) { out.close(); out2.close(); if (pDataArray->trim) { out3.close(); } inFASTA.close(); delete chimera;  return 0;	}
		
		if (chimera->getUnaligned()) { 
			pDataArray->m->mothurOut("Your template sequences are different lengths, please correct."); pDataArray->m->mothurOutEndLine(); 
			out.close(); out2.close(); if (pDataArray->trim) { out3.close(); } inFASTA.close();
			delete chimera;
			return 0; 
		}
		int templateSeqsLength = chimera->getLength();
		
		if (pDataArray->start == 0) { chimera->printHeader(out); }
		
		pDataArray->count = 0;
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
				pDataArray->count++;
			}
			
			delete candidateSeq;
			//report progress
			if((pDataArray->count) % 100 == 0){	pDataArray->m->mothurOutJustToScreen("Processing sequence: " + toString(pDataArray->count) +"\n"); 	}
		}
		//report progress
		if((pDataArray->count) % 100 != 0){	pDataArray->m->mothurOutJustToScreen("Processing sequence: " + toString(pDataArray->count)+"\n"); 		}
		
		pDataArray->numNoParents = chimera->getNumNoParents();
		if (pDataArray->numNoParents == pDataArray->count) { 	pDataArray->m->mothurOut("[WARNING]: megablast returned 0 potential parents for all your sequences. This could be due to formatdb.exe not being setup properly, please check formatdb.log for errors.\n"); }

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

/**************************************************************************************************/

static DWORD WINAPI MySlayerGroupThreadFunction(LPVOID lpParam){ 
	slayerData* pDataArray;
	pDataArray = (slayerData*)lpParam;
	
	try {
        ofstream outCountList;
        if (pDataArray->hasCount && pDataArray->dups) { pDataArray->m->openOutputFile(pDataArray->countlist, outCountList); }

		int totalSeqs = 0;
        pDataArray->end = 0;
		
		for (map<string, map<string, int> >::iterator itFile = pDataArray->fileToPriority.begin(); itFile != pDataArray->fileToPriority.end(); itFile++) {
			
			if (pDataArray->m->control_pressed) {  return 0;  }
			
			int start = time(NULL);
			string thisFastaName = itFile->first;
			map<string, int> thisPriority = itFile->second;
			string thisoutputFileName = pDataArray->m->getRootName(pDataArray->m->getSimpleName(thisFastaName)) + pDataArray->fileGroup[thisFastaName] + "slayer.chimera";
			string thisaccnosFileName = pDataArray->m->getRootName(pDataArray->m->getSimpleName(thisFastaName)) + pDataArray->fileGroup[thisFastaName] + "slayer.accnos";
			string thistrimFastaFileName = pDataArray->m->getRootName(pDataArray->m->getSimpleName(thisFastaName)) + pDataArray->fileGroup[thisFastaName] + "slayer.fasta";
			
			pDataArray->m->mothurOutEndLine(); pDataArray->m->mothurOut("Checking sequences from group: " + pDataArray->fileGroup[thisFastaName] + "."); pDataArray->m->mothurOutEndLine(); 
		
			//int numSeqs = driver(lines[0], thisoutputFileName, thisFastaName, thisaccnosFileName, thistrimFastaFileName, thisPriority);
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			
			ofstream out;
			pDataArray->m->openOutputFile(thisoutputFileName, out);
			
			ofstream out2;
			pDataArray->m->openOutputFile(thisaccnosFileName, out2);
			
			ofstream out3;
			if (pDataArray->trim) {  pDataArray->m->openOutputFile(thistrimFastaFileName, out3); }
			
			ifstream inFASTA;
			pDataArray->m->openInputFile(thisFastaName, inFASTA);
			
			Chimera* chimera;
			chimera = new ChimeraSlayer(thisFastaName, pDataArray->templatefile, pDataArray->trim, thisPriority, pDataArray->search, pDataArray->ksize, pDataArray->match, pDataArray->mismatch, pDataArray->window, pDataArray->divR, pDataArray->minSimilarity, pDataArray->minCoverage, pDataArray->minBS, pDataArray->minSNP, pDataArray->parents, pDataArray->iters, pDataArray->increment, pDataArray->numwanted, pDataArray->realign, pDataArray->blastlocation, pDataArray->threadId);	
			chimera->printHeader(out); 
			
			int numSeqs = 0;
			
			if (pDataArray->m->control_pressed) { out.close(); out2.close(); if (pDataArray->trim) { out3.close(); } inFASTA.close(); delete chimera;  return 0;	}
			
			if (chimera->getUnaligned()) { 
				pDataArray->m->mothurOut("Your template sequences are different lengths, please correct."); pDataArray->m->mothurOutEndLine(); 
				out.close(); out2.close(); if (pDataArray->trim) { out3.close(); } inFASTA.close();
				delete chimera;
				return 0; 
			}
			int templateSeqsLength = chimera->getLength();
			
			bool done = false;
			while (!done) {
				
				if (pDataArray->m->control_pressed) {	out.close(); out2.close(); if (pDataArray->trim) { out3.close(); } inFASTA.close(); delete chimera; return 1;	}
				
				Sequence* candidateSeq = new Sequence(inFASTA);  pDataArray->m->gobble(inFASTA);
				string candidateAligned = candidateSeq->getAligned();
				
				if (candidateSeq->getName() != "") { //incase there is a commented sequence at the end of a file
					if (candidateSeq->getAligned().length() != templateSeqsLength) {  
						pDataArray->m->mothurOut(candidateSeq->getName() + " is not the same length as the template sequences. Skipping."); pDataArray->m->mothurOutEndLine();
					}else{
						//find chimeras
						chimera->getChimeras(candidateSeq);
						
						if (pDataArray->m->control_pressed) {	out.close(); out2.close(); if (pDataArray->trim) { out3.close(); } inFASTA.close(); delete candidateSeq; delete chimera; return 1;	}
						
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
					numSeqs++;
				}
				
				delete candidateSeq;
				
				if (inFASTA.eof()) { break; }
				
				//report progress
				if((numSeqs) % 100 == 0){	pDataArray->m->mothurOutJustToScreen("Processing sequence: " + toString(numSeqs)+"\n"); pDataArray->m->mothurOutEndLine();		}
			}
			//report progress
			if((numSeqs) % 100 != 0){	pDataArray->m->mothurOutJustToScreen("Processing sequence: " + toString(numSeqs)+"\n"); 		}
			
			pDataArray->numNoParents = chimera->getNumNoParents();
			if (pDataArray->numNoParents == numSeqs) { 	pDataArray->m->mothurOut("[WARNING]: megablast returned 0 potential parents for all your sequences. This could be due to formatdb.exe not being setup properly, please check formatdb.log for errors.\n"); }
			
			out.close();
			out2.close();
			if (pDataArray->trim) { out3.close(); }
			inFASTA.close();
			delete chimera;
			pDataArray->end++;
			
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			
            //if we provided a count file with group info and set dereplicate=t, then we want to create a *.pick.count_table
            //This table will zero out group counts for seqs determined to be chimeric by that group.
            if (pDataArray->dups) {
                if (!pDataArray->m->isBlank(thisaccnosFileName)) {
                    ifstream in;
                    pDataArray->m->openInputFile(thisaccnosFileName, in);
                    string name;
                    if (pDataArray->hasCount) {
                        while (!in.eof()) {
                            in >> name; pDataArray->m->gobble(in);
                            outCountList << name << '\t' << pDataArray->fileGroup[thisFastaName] << endl;
                        }
                        in.close();
                    }else {
                        map<string, map<string, string> >::iterator itGroupNameMap = pDataArray->group2NameMap.find(pDataArray->fileGroup[thisFastaName]);
                        if (itGroupNameMap != pDataArray->group2NameMap.end()) {
                            map<string, string> thisnamemap = itGroupNameMap->second;
                            map<string, string>::iterator itN;
                            ofstream out;
                            pDataArray->m->openOutputFile(thisaccnosFileName+".temp", out);
                            while (!in.eof()) {
                                in >> name; pDataArray->m->gobble(in);
                                //pDataArray->m->mothurOut("here = " + name + '\t');
                                itN = thisnamemap.find(name);
                                if (itN != thisnamemap.end()) {
                                    vector<string> tempNames; pDataArray->m->splitAtComma(itN->second, tempNames);
                                    for (int j = 0; j < tempNames.size(); j++) { out << tempNames[j] << endl; }
                                    //pDataArray->m->mothurOut(itN->second + '\n');
                                    
                                }else { pDataArray->m->mothurOut("[ERROR]: parsing cannot find " + name + ".\n"); pDataArray->m->control_pressed = true; }
                            }
                            out.close();
                            in.close();
                            pDataArray->m->renameFile(thisaccnosFileName+".temp", thisaccnosFileName);
                        }else { pDataArray->m->mothurOut("[ERROR]: parsing cannot find " + pDataArray->fileGroup[thisFastaName] + ".\n"); pDataArray->m->control_pressed = true; }
                    }
                    
                }
            }
            
            
			//append files
			pDataArray->m->appendFiles(thisoutputFileName, pDataArray->outputFName); pDataArray->m->mothurRemove(thisoutputFileName); 
			pDataArray->m->appendFiles(thisaccnosFileName, pDataArray->accnos); pDataArray->m->mothurRemove(thisaccnosFileName);
			if (pDataArray->trim) { pDataArray->m->appendFiles(thistrimFastaFileName, pDataArray->fasta); pDataArray->m->mothurRemove(thistrimFastaFileName); }
			pDataArray->m->mothurRemove(thisFastaName);
			
			totalSeqs += numSeqs;
			
			pDataArray->m->mothurOutEndLine(); pDataArray->m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences from group " + pDataArray->fileGroup[thisFastaName] + ".");	pDataArray->m->mothurOutEndLine();
		}
		
		pDataArray->count = totalSeqs;
        if (pDataArray->hasCount && pDataArray->dups) { outCountList.close(); }
		
		return 0;
		
	}
	catch(exception& e) {
		pDataArray->m->errorOut(e, "ChimeraSlayerCommand", "MySlayerGroupThreadFunction");
		exit(1);
	}
} 

#endif

/**************************************************************************************************/


#endif


