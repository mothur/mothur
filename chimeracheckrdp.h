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
#include "database.hpp"

//This class was created using the algorythms described in 
//CHIMERA_CHECK version 2.7 written by Niels Larsen. 

/***********************************************************/

class ChimeraCheckRDP : public Chimera {
	
	public:
		ChimeraCheckRDP(string, string);	
		~ChimeraCheckRDP();
		
		void getChimeras();
		void print(ostream&);
		
		void setCons(string){};
		void setQuantiles(string q) {};
		
		
	private:
		
		vector<linePair*> lines;
		vector<Sequence*> querySeqs;
		Database* templateDB;
		Kmer* kmer;
		
		vector< vector<sim> > IS;  //IS[0] is the vector of IS values for each window for querySeqs[0]
		
		//map of vector of maps- I know its a little convaluted but I am trying to save time 
		//I think that since the window is only sliding 10 bases there is a good probability that the closest seq to each fragment
		//will be the same for several windows so I want to save the vector of maps containing its kmer info rather than regenerating it.
		//So...
		map<string, vector< map<int, int> > > seqKmerInfo;  // outer map - sequence name -> kmer info 
															// kmer info: inner vector of maps - each entry in the vector is a map of the kmers up to that spot in the unaligned seq
															//example:  seqKmerInfo["ecoli"][50] = map containing the kmers found in the first 50 + kmersize characters of ecoli.
															//i chose to store the kmers numbers in a map so you wouldn't have to check for dupilcate entries and could easily find the 
															//kmers 2 seqs had in common.  There may be a better way to do this thats why I am leaving so many comments...
		map<string, vector< map<int, int> > >:: iterator it;
		map<int, int>::iterator it2;
		
		vector<Sequence> closest;		//closest[0] is the closest overall seq to querySeqs[0].
		
		string fastafile, templateFile;
		
		
		vector<sim> findIS(int);
		int calcKmers(map<int, int>, map<int, int>);		
		vector< vector<sim> > createProcessesIS(vector<Sequence*>, vector<linePair*>);
		
};

/***********************************************************/

#endif

