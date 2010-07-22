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
		string leftParent; //keep the name of closest left 
		string rightParent; //keep the name of closest 
		float score;  //preference score
		float closestLeft;  //keep the closest left 
		float closestRight; //keep the closest right 
		int midpoint;
		Preference() { name = ""; leftParent = ""; rightParent = ""; score = 0.0; closestLeft = 10000.0; closestRight = 10000.0; midpoint = 0;  }
		~Preference() {}
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
	string parentAligned;
	float queryToParent;
	float queryToParentLocal;
	float divR;
};
/***********************************************************************/
struct SeqDist {
	Sequence* seq;
	float dist;
	int index;
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
			unsigned long int start;
			unsigned long int end;
			linePair(unsigned long int i, unsigned long int j) : start(i), end(j) {}
			linePair(){}
};

/***********************************************************************/

class Chimera {

	public:
	
		Chimera(){ m = MothurOut::getInstance(); length = 0; unaligned = false; }
		virtual ~Chimera(){	for (int i = 0; i < templateSeqs.size(); i++) { delete templateSeqs[i];  } };
		virtual bool getUnaligned()				{	return unaligned;			}
		virtual int getLength()					{   return length;	}
		virtual vector<Sequence*> readSeqs(string);
		virtual void setMask(string);
		virtual map<int, int> runFilter(Sequence*);
		virtual string createFilter(vector<Sequence*>, float);
		virtual void printHeader(ostream&){};
		virtual int getChimeras(Sequence*){ return 0; }
		virtual int getChimeras(){ return 0; }
		virtual int print(ostream&, ostream&){  return 0; }
		
		#ifdef USE_MPI
		virtual int print(MPI_File&, MPI_File&){  return 0; }
		#endif
		
		
	protected:
		
		vector<Sequence*> templateSeqs;
		bool filter, unaligned; 
		int length; 
		string seqMask, filterString, outputDir, templateFileName; 
		Sequence* getSequence(string);  //find sequence from name	
		MothurOut* m;
};

/***********************************************************************/

#endif

