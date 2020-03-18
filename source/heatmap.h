#ifndef HEATMAP_H
#define HEATMAP_H
/*
 *  heatmap.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/25/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "sharedrabundvectors.hpp"
#include "rabundvector.hpp"
#include "sharedrabundfloatvectors.hpp"
#include "utils.hpp"

/***********************************************************************/
struct binCount {
		int bin;
		int abund;
		binCount(int i, int j) : bin(i), abund(j) {}
};
/***********************************************************************/
struct binCountFloat {
		int bin;
		float abund;
		binCountFloat(int i, float j) : bin(i), abund(j) {}
};

/***********************************************************************/
//sorts highest abund to lowest
inline bool comparebinCounts(binCount left, binCount right){
	return (left.abund > right.abund);	
}
/***********************************************************************/
//sorts highest abund to lowest
inline bool comparebinFloatCounts(binCountFloat left, binCountFloat right){
	return (left.abund > right.abund);	
}
/***********************************************************************/

class HeatMap {
	
	public:
		HeatMap(string, string, int, int, string, string);
		~HeatMap(){};
	
		string getPic(RAbundVector*);
		string getPic(SharedRAbundVectors*&);
		string getPic(SharedRAbundFloatVectors*&);

	private:
		vector<string> sortSharedVectors(vector<SharedRAbundVector*>, vector<string>);
		vector<string> sortSharedVectors(vector<SharedRAbundFloatVector*>, vector<string>);
		int sortRabund(RAbundVector*);
		void printLegend(int, float);

		string format, sorted, groupComb, scaler, outputDir, inputfile;
		ofstream outsvg;
		MothurOut* m;
        Utils util;
		int numOTU, fontSize;
		
		map<int, int> orderTopGroup(vector<SharedRAbundVector*>);
		map<int, int> orderTopOtu(vector<SharedRAbundVector*>);
		map<int, int> orderShared(vector<SharedRAbundVector*>);
		map<int, int> orderTopGroup(vector<SharedRAbundFloatVector*>);
        map<int, int> orderTopOtu(vector<SharedRAbundFloatVector*>);
		map<int, int> orderShared(vector<SharedRAbundFloatVector*>);

			
};

/***********************************************************************/

#endif




