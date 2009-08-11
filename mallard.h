#ifndef MALLARD_H
#define MALLARD_H

/*
 *  mallard.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 8/11/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


#include "chimera.h"
#include "decalc.h"

//This class was created using the algorythms described in the 
// "New Screening Software Shows that Most Recent Large 16SrRNA Gene Clone Libraries Contain Chimeras" paper 
//by Kevin E. Ashelford 1, Nadia A. Chuzhanova 2, John C. Fry 1, Antonia J. Jones 3 and Andrew J. Weightman 1.

/***********************************************************/

class Mallard : public Chimera {
	
	public:
		Mallard(string);	
		~Mallard();
		
		void getChimeras();
		void print(ostream&);
		
		void setCons(string c)	{};
		void setQuantiles(string q) {};
		
	private:
	
		struct linePair {
			int start;
			int end;
			linePair(int i, int j) : start(i), end(j) {}
			linePair(){}
		};

		DeCalculator* decalc;
		int iters;
		string fastafile;
		
		vector<linePair*> lines;
		vector<Sequence*> querySeqs;
		vector<int> windowSizes;			//windowSizes[0] = window size of querySeqs[0]
		vector<int> marked;
		vector<float> highestDE;
	
		vector<float> probabilityProfile;
		vector< vector<quanMember> > quantilesMembers;  //quantiles[0] is the vector of deviations with ceiling score of 1, quantiles[1] is the vector of deviations with ceiling score of 2...
		
		void createProcessesQuan();
		
};

/***********************************************************/

#endif

