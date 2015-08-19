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
struct data_struct { 
	float divr_qla_qrb;
	float divr_qlb_qra;
	float qla_qrb;
	float qlb_qra;
	float qla;
	float qrb;
	float ab; 
	float qa;
	float qb; 
	float lab; 
	float rab; 
	float qra; 
	float qlb; 
	int winLStart;
	int winLEnd; 
	int winRStart; 
	int winREnd; 
	Sequence querySeq; 
	Sequence parentA;
	Sequence parentB;
	float bsa;
	float bsb;
	float bsMax;
	float chimeraMax;
	
};
/***********************************************************************/
struct data_results {
	vector<data_struct> results;
	string flag;
	Sequence trimQuery;
	//results malignerResults;
	
	data_results(vector<data_struct> d, string f, map<int, int> s, Sequence t) : results(d), flag(f), trimQuery(t) {}
	data_results() {}
};
/***********************************************************************/
//sorts lowest to highest first by bsMax, then if tie by chimeraMax
inline bool compareDataStruct(data_struct left, data_struct right){
	if (left.bsMax < right.bsMax) { return true; }
	else if (left.bsMax == right.bsMax) {
		return (left.chimeraMax < right.chimeraMax);
	}else { return false;	}
} 
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
//	int mismatches;
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
/***********************************************************************/
struct SeqCompare {
	Sequence seq;
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
//sorts lowest to highest
inline bool compareSeqCompare(SeqCompare left, SeqCompare right){
	return (left.dist < right.dist);	
} 
//********************************************************************************************************************
struct sim {
		string leftParent;
		string rightParent; 
		float score;  
		int midpoint;
};
/***********************************************************************/

class Chimera {

	public:
	
		Chimera(){ m = MothurOut::getInstance(); length = 0; unaligned = false;  byGroup = false; }
		virtual ~Chimera(){	for (int i = 0; i < templateSeqs.size(); i++) { delete templateSeqs[i];  } for (int i = 0; i < filteredTemplateSeqs.size(); i++) { delete filteredTemplateSeqs[i];  } };
		virtual bool getUnaligned()				{	return unaligned;			}
		virtual int getLength()					{   return length;	}
		virtual vector<Sequence*> readSeqs(string);
		virtual void setMask(string);
		virtual map<int, int> runFilter(Sequence*);
		virtual string createFilter(vector<Sequence*>, float);
		virtual void printHeader(ostream&){};
		virtual int getChimeras(Sequence*){ return 0; }
		virtual int getChimeras(){ return 0; }
		virtual Sequence print(ostream&, ostream&){  Sequence temp; return temp; }
		virtual Sequence print(ostream&, ostream&, data_results, data_results) { Sequence temp; return temp; }
		virtual int print(ostream&, ostream&, string){  return 0; }
		virtual int getNumNoParents(){  return 0; }
		virtual data_results getResults() { data_results results; return results; }
				
	protected:
		
		vector<Sequence*> templateSeqs;
		vector<Sequence*> filteredTemplateSeqs;
		bool filter, unaligned, byGroup; 
		int length; 
		string seqMask, filterString, outputDir, templateFileName; 
		Sequence* getSequence(string);  //find sequence from name	
		MothurOut* m;
};

/***********************************************************************/

#endif

