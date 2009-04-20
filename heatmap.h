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

#include "ordervector.hpp"
#include "rabundvector.hpp"
#include "sharedrabundvector.h"
#include "sharedordervector.h"
#include "datavector.hpp"
#include "globaldata.hpp"
#include "sharedutilities.h"

/***********************************************************************/

class HeatMap {
	
	public:
		HeatMap();
		~HeatMap(){ delete util; };
	
		void getPic(OrderVector*);
		void getPic(SharedOrderVector*);

	private:
		void sortSharedVectors(vector<SharedRAbundVector*>& );
		void printLegend(int, float);

		GlobalData* globaldata;
		SharedUtil* util;
		string format, sorted, groupComb, scaler;
		ofstream outsvg;
			
};

/***********************************************************************/

#endif



