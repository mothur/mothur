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
#include "sequence.hpp"
#include "dist.h"

struct Preference {
		string name;
		vector<string> leftParent; //keep the name of closest left associated with the two scores
		vector<string> rightParent; //keep the name of closest right associated with the two scores
		vector<float> score;  //so you can keep last score and calc this score and keep whichever is bigger.
		vector<float> closestLeft;  //keep the closest left associated with the two scores
		vector<float> closestRight; //keep the closest right associated with the two scores
		int midpoint;

};


/***********************************************************/

class Bellerophon : public Chimera {
	
	public:
		Bellerophon(string);	
		~Bellerophon() {};
		
		void getChimeras();
		void print(ostream&);
		
		void setCons(string){};
		void setQuantiles(string) {};
		
		
	private:
		Dist* distCalculator;
		FilterSeqsCommand* filterSeqs;
		vector<Sequence> seqs;
		vector<Preference> pref;
		string fastafile;
		int iters;
		
		void generatePreferences(vector<SeqMap>, vector<SeqMap>, int);
		int createSparseMatrix(int, int, SparseMatrix*, vector<Sequence>);
};

/***********************************************************/

#endif

