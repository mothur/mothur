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

/***********************************************************************/
struct Preference {
		string name;
		vector<string> leftParent; //keep the name of closest left associated with the two scores
		vector<string> rightParent; //keep the name of closest right associated with the two scores
		vector<float> score;  //so you can keep last score and calc this score and keep whichever is bigger.
		vector<float> closestLeft;  //keep the closest left associated with the two scores
		vector<float> closestRight; //keep the closest right associated with the two scores
		int midpoint;

};
/***********************************************************************/
struct score_struct {
	int prev;
	int score;
	int row;
	int col;
};
/***********************************************************************/
struct trace_struct {
	int col;
	int oldCol;
	int row;
};
/***********************************************************************/
struct results {
	int regionStart;
	int regionEnd;
	int nastRegionStart;
	int nastRegionEnd;
	string parent;
	float queryToParent;
	float queryToParentLocal;
	float divR;
};
/***********************************************************************/
struct SeqDist {
	Sequence* seq;
	float dist;
};
//********************************************************************************************************************
//sorts lowest to highest
inline bool compareRegionStart(results left, results right){
	return (left.nastRegionStart < right.nastRegionStart);	
} 
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
		Chimera(string, bool);
		Chimera(string, string);
		virtual ~Chimera(){};
		virtual void setFilter(bool f)			{	filter = f;	 		}
		virtual void setCorrection(bool c)		{	correction = c;		}
		virtual void setProcessors(int p)		{	processors = p;		}
		virtual void setWindow(int w)			{	window = w;			}
		virtual void setIncrement(int i)		{	increment = i;		}
		virtual void setNumWanted(int n)		{	numWanted = n;		}
		virtual void setKmerSize(int k)			{	kmerSize = k;		}
		virtual void setSVG(int s)				{	svg = s;			}
		virtual void setName(string n)			{	name = n;			}
		virtual void setMatch(int m)			{	match = m;			}
		virtual void setMisMatch(int m)			{	misMatch = m;		}
		virtual void setDivR(float d)			{	divR = d;			}
		virtual void setParents(int p)			{	parents = p;		}
		virtual void setMinSim(int s)			{	minSim = s;			}
		virtual void setMinCoverage(int c)		{	minCov = c;			}
		virtual void setMinBS(int b)			{	minBS = b;			}
		virtual void setMinSNP(int s)			{	minSNP = s;			}
		virtual void setIters(int i)			{	iters = i;			}
		virtual void setTemplateSeqs(vector<Sequence*> t)	{	templateSeqs = t;	}
		virtual bool getUnaligned()				{	return unaligned;			}
		virtual void setTemplateFile(string t)	{   templateFileName = t;	}
		
		virtual void setCons(string){};
		virtual void setQuantiles(string){};
		virtual void doPrep(){};
		virtual vector<Sequence*> readSeqs(string);
		virtual vector< vector<float> > readQuantiles();
		virtual void setMask(string);
		virtual void runFilter(Sequence*);
		virtual void createFilter(vector<Sequence*>);
		
		virtual void printHeader(ostream&){};
		virtual int getChimeras(Sequence*){ return 0; }
		virtual int getChimeras(){ return 0; }
		virtual void print(ostream&){};	
		
		
	protected:
		
		vector<Sequence*> templateSeqs;
		bool filter, correction, svg, unaligned;
		int processors, window, increment, numWanted, kmerSize, match, misMatch, minSim, minCov, minBS, minSNP, parents, iters;
		float divR;
		string seqMask, quanfile, filterString, name, outputDir, templateFileName;
		Sequence* getSequence(string);  //find sequence from name	
};

/***********************************************************************/

#endif

