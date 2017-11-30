#ifndef SEQNOISE
#define SEQNOISE



/*
 *  mySeqNoise.h
 *  
 *
 *  Created by Pat Schloss on 8/31/11.
 *  Copyright 2011 Patrick D. Schloss. All rights reserved.
 *
 */

/*****************************************************************************************************************************/
/*****************************************************************************************************************************/
/* NOTE: Order matters in this class.  If you are going to use it, make sure your files have the sequences in the same order. */
/*****************************************************************************************************************************/
/*****************************************************************************************************************************/


#include "mothurout.h"
#include "utils.hpp"
/**************************************************************************************************/

struct freqData {
	
	freqData(int i, int freq) : frequency(freq), index(i){	}
	
	bool operator<( freqData const& rhs ) const {
		return frequency < rhs.frequency; 
	}
	
	int frequency;
	int index;
	
};
/**************************************************************************************************/

class seqNoise {
public:
	seqNoise() { m = MothurOut::getInstance();  }
	~seqNoise(){}
	
	int getSequenceData(string, vector<string>&);
	int addSeq(string, vector<string>&); 
	int getRedundantNames(string, vector<string>&, vector<string>&, vector<int>&);
	int addRedundantName(string, string, vector<string>&, vector<string>&, vector<int>&);
    int getDistanceData(string, vector<double>&);
	int getListData(string, double, vector<int>&, vector<int>&, vector<vector<int> >&);
	int updateOTUCountData(vector<int>, vector<vector<int> >, vector<vector<int> >, vector<int>&, vector<int>&, vector<int>&);
	double calcNewWeights(vector<double>&,vector<int>,vector<int>,vector<int>,vector<int>,vector<int>,vector<double>);
	int calcCentroids(vector<int>,vector<int>,vector<int>&,vector<int>&,vector<int>,vector<double>,vector<int>,vector<int>,vector<double>);
	int checkCentroids(vector<double>&, vector<int>);
	int setUpOTUData(vector<int>&, vector<double>&, vector<int>, vector<double>, vector<int>, vector<int>, vector<int>);
	int finishOTUData(vector<int>, vector<int>&, vector<int>&, vector<int>&, vector<int>&, vector<vector<int> >&, vector<vector<int> >&, vector<double>&);
	int writeOutput(string, string, string, vector<int>, vector<int>, vector<int>, vector<string>, vector<string>, vector<string>, vector<int>, vector<double>&);


private:
	MothurOut* m;
    Utils util;
	
	int getLastMatch(char, vector<vector<char> >&, int, int, vector<int>&, vector<int>&);
	int countDiffs(vector<int>, vector<int>);
	vector<int> convertSeq(string);
	string degapSeq(string);
	
};

/**************************************************************************************************/
#endif

