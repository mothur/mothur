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


#include "chimera.h"
#include "maligner.h"
#include "slayer.h"

//***********************************************************************/
//This class was modeled after the chimeraSlayer written by the Broad Institute
/***********************************************************************/


class ChimeraSlayer : public Chimera {
	
	public:
		ChimeraSlayer(string, string, string, int, int, int, int, float, int, int, int, int, int, int, int, int, bool);
		~ChimeraSlayer();
		
		int getChimeras(Sequence*);
		int print(ostream&, ostream&);
		void printHeader(ostream&);
		int doPrep();
		
		#ifdef USE_MPI
		int print(MPI_File&, MPI_File&);
		#endif
		
	private:
		Sequence* querySeq;
		DeCalculator* decalc;
		Maligner* maligner;
		Slayer* slayer;
		map<int, int>  spotMap;
		Database* databaseRight;
		Database* databaseLeft;
		
		vector<data_struct>  chimeraResults;
		string chimeraFlags, searchMethod, fastafile;
		bool realign;
		int window, numWanted, kmerSize, match, misMatch, minSim, minCov, minBS, minSNP, parents, iters, increment;
		float divR;
	
		void printBlock(data_struct, ostream&);
		string getBlock(data_struct);
		
};

/************************************************************************/

#endif


