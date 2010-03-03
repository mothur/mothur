#ifndef BELLEROPHON_H
#define BELLEROPHON_H

/*
 *  bellerophon.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


#include "chimera.h"
#include "filterseqscommand.h"
#include "sparsematrix.hpp"
#include "sequence.hpp"
#include "dist.h"

typedef list<PCell>::iterator MatData;
typedef map<int, float> SeqMap;  //maps sequence to all distance for that seqeunce

/***********************************************************/

class Bellerophon : public Chimera {
	
	public:
		Bellerophon(string, string);	
		~Bellerophon() {};
		
		int getChimeras();
		int print(ostream&, ostream&);
		
	private:
		Dist* distCalculator;
		FilterSeqsCommand* filterSeqs;
		vector<Sequence*> seqs;
		vector<Preference> pref;
		string fastafile;
		int iters;
		
		int generatePreferences(vector<SeqMap>, vector<SeqMap>, int);
		int createSparseMatrix(int, int, SparseMatrix*, vector<Sequence>);
};

/***********************************************************/

#endif

