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

/***********************************************************************/
//This class was modeled after the chimeraSlayer written by the Broad Institute
/***********************************************************************/


class ChimeraSlayer : public Chimera {
	
	public:
		ChimeraSlayer(string, bool, string);	
		~ChimeraSlayer();
		
		int getChimeras(Sequence*);
		void print(ostream&, ostream&);
		void printHeader(ostream&);
		void doPrep();
		
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
	
		void printBlock(data_struct, ostream&);
		
};

/************************************************************************/

#endif


