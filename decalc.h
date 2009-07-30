#ifndef DECALC_H
#define DECALC_H
/*
 *  decalc.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/22/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "mothur.h"
#include "sequence.hpp"

/***********************************************************************/

//This class was created using the algorythms described in the 
// "At Least 1 in 20 16S rRNA Sequence Records Currently Held in the Public Repositories is Estimated To Contain Substantial Anomalies" paper 
//by Kevin E. Ashelford 1, Nadia A. Chuzhanova 3, John C. Fry 1, Antonia J. Jones 2 and Andrew J. Weightman 1.

/***********************************************************************/

class DeCalculator {

	public:
		
		DeCalculator() {};
		~DeCalculator() {};
		
		set<int> getPos() {  return h;  }
		void setMask(string m); 
		void runMask(Sequence*);
		void trimSeqs(Sequence*, Sequence*, map<int, int>&);
		void removeObviousOutliers(vector< vector<float> >&);
		vector<float> calcFreq(vector<Sequence*>, string);
		vector<int> findWindows(Sequence*, int, int, int&, int);
		vector<float> calcObserved(Sequence*, Sequence*, vector<int>, int);
		vector<float>  calcExpected(vector<float>, float);
		vector<float>  findQav(vector<int>, int, vector<float>);  
		float calcDE(vector<float>, vector<float>);
		float calcDist(Sequence*, Sequence*, int, int);
		float getCoef(vector<float>, vector<float>);
		vector< vector<float> > getQuantiles(vector<Sequence*>, vector<int>, int, vector<float>, int, int, int);
		
	private:
		float findAverage(vector<float> myVector);
		string seqMask;
		set<int> h;

};

/***********************************************************************/

#endif
