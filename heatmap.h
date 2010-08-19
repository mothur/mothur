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
#include "sharedrabundvector.h"
#include "datavector.hpp"
#include "globaldata.hpp"

/***********************************************************************/
struct binCount {
		int bin;
		int abund;
		binCount(int i, int j) : bin(i), abund(j) {}
};
/***********************************************************************/
//sorts highest abund to lowest
inline bool comparebinCounts(binCount left, binCount right){
	return (left.abund > right.abund);	
}
/***********************************************************************/

class HeatMap {
	
	public:
		HeatMap(string, string, int, int, string);
		~HeatMap(){};
	
		string getPic(RAbundVector*);
		string getPic(vector<SharedRAbundVector*>);

	private:
		int sortSharedVectors(vector<SharedRAbundVector*>& );
		void printLegend(int, float);

		GlobalData* globaldata;
		string format, sorted, groupComb, scaler, outputDir;
		ofstream outsvg;
		MothurOut* m;
		int numOTU, fontSize;
		
		map<int, int> orderTopGroup(vector<SharedRAbundVector*>&);
		map<int, int> orderTopOtu(vector<SharedRAbundVector*>&);
		map<int, int> orderShared(vector<SharedRAbundVector*>&);
			
};

/***********************************************************************/

#endif




