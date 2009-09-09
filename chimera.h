#ifndef CHIMERA_H
#define CHIMERA_H

/*
 *  chimera.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


#include "mothur.h"
#include "sequence.hpp"


struct Preference {
		string name;
		vector<string> leftParent; //keep the name of closest left associated with the two scores
		vector<string> rightParent; //keep the name of closest right associated with the two scores
		vector<float> score;  //so you can keep last score and calc this score and keep whichever is bigger.
		vector<float> closestLeft;  //keep the closest left associated with the two scores
		vector<float> closestRight; //keep the closest right associated with the two scores
		int midpoint;

};

struct SeqDist {
	Sequence* seq;
	float dist;
};

//********************************************************************************************************************
//sorts lowest to highest
inline bool compareSeqDist(SeqDist left, SeqDist right){
	return (left.dist < right.dist);	
} 
//********************************************************************************************************************

struct sim {
		string leftParent;
		string rightParent; 
		float score;  
		int midpoint;
};

struct linePair {
			int start;
			int end;
			linePair(int i, int j) : start(i), end(j) {}
			linePair(){}
};

/***********************************************************************/

class Chimera {

	public:
	
		Chimera(){};
		Chimera(string);
		Chimera(string, string);
		virtual ~Chimera(){};
		virtual void setFilter(bool f)			{	filter = f;	 		}
		virtual void setCorrection(bool c)		{	correction = c;		}
		virtual void setProcessors(int p)		{	processors = p;		}
		virtual void setWindow(int w)			{	window = w;			}
		virtual void setIncrement(int i)		{	increment = i;		}
		virtual void setNumWanted(int n)		{	numWanted = n;		}
		virtual void setKmerSize(int k)			{	kmerSize = k;		}
		
		virtual void setCons(string){};
		virtual void setQuantiles(string){};
		virtual vector<Sequence*> readSeqs(string);
		virtual vector< vector<float> > readQuantiles();
		virtual void setMask(string);
		virtual void runFilter(vector<Sequence*>);
		virtual void createFilter(vector<Sequence*>);
		
		
		//pure functions
		virtual void getChimeras() = 0;	
		virtual void print(ostream&) = 0;	
		
	protected:
		
		bool filter, correction;
		int processors, window, increment, numWanted, kmerSize;
		string seqMask, quanfile, filterString;
			

};

/***********************************************************************/

#endif

