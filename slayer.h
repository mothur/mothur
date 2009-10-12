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
struct data_struct { //this is crazy big...but follow original.
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
//sorts lowest to highest first by bsMax, then if tie by chimeraMax
inline bool compareDataStruct(data_struct left, data_struct right){
	if (left.bsMax < right.bsMax) { return true; }
	else if (left.bsMax == right.bsMax) {
		return (left.chimeraMax < right.chimeraMax);
	}else { return false;	}
} 
/***********************************************************************/
struct snps { 
	char queryChar;
	char parentAChar;
	char parentBChar;
};

/***********************************************************************/


class Slayer {

	public:
		
		Slayer(int, int, int, float, int);
		~Slayer() {};
		
		string getResults(Sequence*, vector<Sequence*>);
		vector<data_struct> getOutput()  {	return outputResults;			}
		
				
	private:
		
		int windowSize, windowStep, parentFragmentThreshold, iters;
		float divRThreshold; 
		vector<data_struct>  outputResults;
		
		map<int, int> verticalFilter(vector<Sequence*>);
		float computePercentID(string, string, int, int);
		
		vector<data_struct> runBellerophon(Sequence*, Sequence*, Sequence*);
		vector<snps> getSNPS(string, string, string, int, int);
		void bootstrapSNPS(vector<snps>, vector<snps>, float&, float&);
		float snpQA(vector<snps>);
		float snpQB(vector<snps>);
		float snpAB(vector<snps>);
				
};

/***********************************************************************/

#endif


