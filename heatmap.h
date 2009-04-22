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

using namespace std;

#include "rabundvector.hpp"
#include "sharedrabundvector.h"
#include "datavector.hpp"
#include "globaldata.hpp"

/***********************************************************************/

class HeatMap {
	
	public:
		HeatMap();
		~HeatMap(){};
	
		void getPic(RAbundVector*);
		void getPic(vector<SharedRAbundVector*>);

	private:
		void sortSharedVectors(vector<SharedRAbundVector*>& );
		void printLegend(int, float);

		GlobalData* globaldata;
		string format, sorted, groupComb, scaler;
		ofstream outsvg;
			
};

/***********************************************************************/

#endif




