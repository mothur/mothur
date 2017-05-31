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

#include "rabundvector.hpp"
#include "rabundfloatvector.hpp"
#include "datavector.hpp"


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
		string getPic(vector<RAbundVector*>, vector<string>);
		string getPic(vector<RAbundFloatVector*>, vector<string>);

	private:
		vector<string> sortSharedVectors(vector<RAbundVector*>& );
		vector<string> sortSharedVectors(vector<RAbundFloatVector*>& );
		int sortRabund(RAbundVector*&);
		void printLegend(int, float);

		string format, sorted, groupComb, scaler, outputDir, inputfile;
		ofstream outsvg;
		MothurOut* m;
		int numOTU, fontSize;
		
		map<int, int> orderTopGroup(vector<RAbundVector*>&);
		map<int, int> orderTopOtu(vector<RAbundVector*>&);
		map<int, int> orderShared(vector<RAbundVector*>&);
		map<int, int> orderTopGroup(vector<RAbundFloatVector*>&);
		map<int, int> orderTopOtu(vector<RAbundFloatVector*>&);
		map<int, int> orderShared(vector<RAbundFloatVector*>&);

			
};

/***********************************************************************/

#endif




