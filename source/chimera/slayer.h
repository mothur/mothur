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
#include "mothurchimera.h"

/***********************************************************************/
//This class was modeled after the chimeraSlayer written by the Broad Institute
/***********************************************************************/
struct snps { 
	char queryChar;
	char parentAChar;
	char parentBChar;
};

/***********************************************************************/


class Slayer {

	public:
		
		Slayer(int, int, int, float, int, int, int);
		~Slayer() {};
		
		string getResults(Sequence, vector<Sequence>);
		vector<data_struct> getOutput()  {	return outputResults;			}
		
				
	private:
		
		int windowSize, windowStep, parentFragmentThreshold, iters, percentSNPSample, minBS;
		float divRThreshold; 
		vector<data_struct>  outputResults;
		vector< map<int, int> > baseSpots;
		Sequence myQuery;
		
		map<int, int> verticalFilter(Sequence&, Sequence&, Sequence&);
		float computePercentID(string, string, int, int);
		
		vector<data_struct> runBellerophon(Sequence, Sequence, Sequence, map<int, int>&);
		vector<snps> getSNPS(string, string, string, int, int);
		int bootstrapSNPS(vector<snps>, vector<snps>, float&, float&, int);
		float snpQA(vector<snps>);
		float snpQB(vector<snps>);
		float snpAB(vector<snps>);
		MothurOut* m;
        Utils util;
				
};

/***********************************************************************/

#endif


