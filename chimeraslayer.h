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
		ChimeraSlayer(string, string);	
		~ChimeraSlayer();
		
		int getChimeras();
		void print(ostream&);
		
		void setCons(string){};
		void setQuantiles(string q) {};
		
		
	private:
		DeCalculator* decalc;
		Maligner* maligner;
		Slayer* slayer;
		vector<linePair*> lines;
		vector<Sequence*> querySeqs;
		vector<Sequence*> templateSeqs;
		vector< map<int, int> > spotMap;
		
		vector< vector<data_struct> > chimeraResults;
		vector<string> chimeraFlags;
				
		string fastafile, templateFile;
		
		Sequence* getSequence(string);  //find sequence from name
		void printBlock(data_struct, ostream&, int i);
};

/************************************************************************/

#endif


