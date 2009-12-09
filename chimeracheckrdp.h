#ifndef CHIMERACHECK_H
#define CHIMERACHECK_H

/*
 *  chimeracheckrdp.h
 *  Mothur
 *
 *  Created by westcott on 9/8/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */


#include "chimera.h"
#include "kmer.hpp"
#include "kmerdb.hpp"
#include "alignmentdb.h"

/***********************************************************/
//This class was created using the algorythms described in 
//CHIMERA_CHECK version 2.7 written by Niels Larsen. 

/***********************************************************/

class ChimeraCheckRDP : public Chimera {
	
	public:
		ChimeraCheckRDP(string, string);	
		~ChimeraCheckRDP();
		
		int getChimeras();
		void print(ostream&);
		
		void setCons(string){};
		void setQuantiles(string q) {};
		
		
	private:
		
		vector<linePair*> lines;
		vector<Sequence*> querySeqs;
		AlignmentDB* templateDB;
		Kmer* kmer;
		
		vector< vector<sim> > IS;  //IS[0] is the vector of IS values for each window for querySeqs[0]
		float chimeraCutoff;
		
		
		//map<string, vector< map<int, int> > >:: iterator it;
		//map<int, int>::iterator it2;
		
		vector<Sequence> closest;		//closest[0] is the closest overall seq to querySeqs[0].
		
		string fastafile, templateFile;
		map<string, string> names;
		
		vector<sim> findIS(int);
		int calcKmers(map<int, int>, map<int, int>);
		void getCutoff();
		void makeSVGpic(vector<sim>, int);
		void readName(string);
				
		void createProcessesIS();
		
};

/***********************************************************/

#endif

