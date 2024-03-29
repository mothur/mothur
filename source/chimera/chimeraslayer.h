#ifndef CHIMERASLAYER_H
#define CHIMERASLAYER_H

/*
 *  chimeraslayer.h
 *  Mothur
 *
 *  Created by westcott on 9/25/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */


#include "mothurchimera.h"
#include "maligner.h"
#include "slayer.h"



//***********************************************************************/
//This class was modeled after the chimeraSlayer written by the Broad Institute
/***********************************************************************/

class ChimeraSlayer : public MothurChimera {
	
	public:
		ChimeraSlayer(string, string, bool, int, int, int, int, float, int, int, int, int, int, int, int, int, bool, int);
		ChimeraSlayer(string, string, bool, map<string, int>&,  int, int, int, int, float, int, int, int, int, int, int, int, int, bool, int);
		ChimeraSlayer(string, string, bool, map<string, int>&,  int, int, int, int, float, int, int, int, int, int, int, int, int, bool, int, bool);

		~ChimeraSlayer();
		
		int getChimeras(Sequence*);
		Sequence print(ostream&, ostream&);
		Sequence print(ostream&, ostream&, data_results, data_results);
		void printHeader(ostream&);
		int doPrep();
		int getNumNoParents() { return numNoParents; }
		data_results getResults() { return printResults; }
		
	private:
		Sequence querySeq;
		Sequence trimQuery;
		DeCalculator decalc;
        SearchDatabase* databaseRight;
        SearchDatabase* databaseLeft;
		map<string, int> priority; //for template=self, seqname, seqAligned, abundance
		set<string> chimericSeqs; //for template=self, so we don't add chimeric sequences to the userTemplate set
		int numNoParents, threadID;
	
		vector<data_struct>  chimeraResults;
		data_results printResults;
		string chimeraFlags, fastafile;
		bool realign, trimChimera;
		int window, numWanted, kmerSize, match, misMatch, minSim, minCov, minBS, minSNP, parents, iters, increment;
		float divR;
	
		void printBlock(data_struct, string, ostream&);
		void printBlock(data_results, data_results, bool, bool, string, ostream&);
		string getBlock(data_struct, string);
		string getBlock(data_results, data_results, bool, bool, string);
		//int readNameFile(string);
		vector<Sequence*> getTemplate(Sequence, vector<Sequence*>&);
		vector<Sequence> getRefSeqs(Sequence, vector<Sequence*>&, vector<Sequence*>&);
		vector<Sequence> getKmerSeqs(Sequence, vector<Sequence*>&, int);
		
};

/************************************************************************/

#endif


