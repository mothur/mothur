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
//This class was created using the algorithms described in 
//CHIMERA_CHECK version 2.7 written by Niels Larsen. 

/***********************************************************/

class ChimeraCheckRDP : public Chimera {
	
	public:
		ChimeraCheckRDP(string, string, string, bool, int, int, string); //fasta, template, name, svg, increment, ksize, outputDir	
		~ChimeraCheckRDP();
		
		int getChimeras(Sequence*);
		Sequence print(ostream&, ostream&);
		
	private:
		
		Sequence* querySeq;
		AlignmentDB* templateDB;
		Kmer* kmer;
		Sequence closest;		//closest is the closest overall seq to query

		vector<sim>  IS;  //IS is the vector of IS values for each window for query
		string fastafile;
		map<string, string> names;
		string name;
		bool svg;
		int kmerSize, increment;
		
		vector<sim> findIS();
		int calcKmers(map<int, int>, map<int, int>);
		void makeSVGpic(vector<sim>);
		void readName(string);
};
/***********************************************************/

#endif

