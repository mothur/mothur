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
		ChimeraSlayer(string, bool);	
		~ChimeraSlayer();
		
		int getChimeras(Sequence*);
		void print(ostream&);
		void printHeader(ostream&);
		
	private:
		Sequence* querySeq;
		DeCalculator* decalc;
		Maligner* maligner;
		Slayer* slayer;
		map<int, int>  spotMap;
		
		vector<data_struct>  chimeraResults;
		string chimeraFlags, searchMethod;
		string fastafile;
		bool realign;
	
		void printBlock(data_struct, ostream&);
};

/************************************************************************/

#endif


