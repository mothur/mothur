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

//this structure is necessary to determine the sequence that contributed to the outliers when we remove them
//this way we can remove all scores that are contributed by outlier sequences.
struct quanMember {
	float score;
	int member1;
	int member2;
	quanMember (float s, int m, int n) : score(s), member1(m),  member2(n) {}
	quanMember() {}
	
};
		
//********************************************************************************************************************
class DeCalculator {

	public:
		
		DeCalculator() { m = MothurOut::getInstance(); }
		~DeCalculator() {};
		
		vector<Sequence*> findClosest(Sequence*, vector<Sequence*>&, vector<Sequence*>&, int&);  //takes querySeq, a reference db, filteredRefDB, numWanted 
		Sequence* findClosest(Sequence*, vector<Sequence*>);
		set<int> getPos() {  return h;  }
		void setMask(string); 
		void setAlignmentLength(int l) {  alignLength = l;  }
		void runMask(Sequence*);
		void trimSeqs(Sequence*, Sequence*, map<int, int>&);
		map<int, int> trimSeqs(Sequence*, vector<Sequence*>);
		void removeObviousOutliers(vector< vector<float> >&, int);
		vector<float> calcFreq(vector<Sequence*>, string);
		vector<int> findWindows(Sequence*, int, int, int&, int);
		vector<float> calcObserved(Sequence*, Sequence*, vector<int>, int);
		vector<float>  calcExpected(vector<float>, float);
		vector<float>  findQav(vector<int>, int, vector<float>);  
		float calcDE(vector<float>, vector<float>);
		float calcDist(Sequence*, Sequence*, int, int);
		float getCoef(vector<float>, vector<float>);
		vector< vector<float> > getQuantiles(vector<Sequence*>, vector<int>, int, vector<float>, int, int, int);
		
		vector<int> returnObviousOutliers(vector< vector<quanMember> >, int);
		
		map<int, int> getMaskMap() { return maskMap; }
		
	private:
		//vector<quanMember> sortContrib(map<quanMember*, float>);  //used by mallard
		float findAverage(vector<float>);
		//int findLargestContrib(vector<int>);
		//void removeContrib(int, vector<quanMember>&);
		string seqMask;
		set<int> h;
		int alignLength;
		map<int, int> maskMap;
		MothurOut* m;

};

/***********************************************************************/

#endif
