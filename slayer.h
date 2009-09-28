#ifndef SLAYER_H
#define SLAYER_H
/*
 *  slayer.h
 *  Mothur
 *
 *  Created by westcott on 9/25/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

 
#include "sequence.hpp"

/***********************************************************************/
//This class was modeled after the chimeraSlayer written by the Broad Institute
/***********************************************************************/
struct data_struct { //not right needs work...
	int regionStart;
	int regionEnd;
	string parent;
	float queryToParent;
	float queryToParentLocal;
	float divR;
};
/***********************************************************************/


class Slayer {

	public:
		
		Slayer(int, int, int, float);
		~Slayer() {};
		
		void getResults(Sequence*, vector<Sequence*>);
		//float getPercentID() {	return percentIdenticalQueryChimera;	}
		//vector<results> getOutput()  {	return outputResults;			}
		
				
	private:
		
		int windowSize, windowStep, parentFragmentThreshold;
		float divRThreshold; 
		
		void verticalFilter(vector<Sequence*>);
		float computePercentID(string, string, int, int);
		
		vector<data_struct> runBellerophon(Sequence*, Sequence*, Sequence*);
		
				
};

/***********************************************************************/

#endif


