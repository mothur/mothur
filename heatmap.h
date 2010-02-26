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

class HeatMap {
	
	public:
		HeatMap(string, string, string);
		~HeatMap(){};
	
		string getPic(RAbundVector*);
		string getPic(vector<SharedRAbundVector*>);

	private:
		void sortSharedVectors(vector<SharedRAbundVector*>& );
		void printLegend(int, float);

		GlobalData* globaldata;
		string format, sorted, groupComb, scaler, outputDir;
		ofstream outsvg;
		MothurOut* m;
			
};

/***********************************************************************/

#endif




