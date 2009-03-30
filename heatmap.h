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

/***********************************************************************/

class HeatMap {
	
	public:
		HeatMap();
		~HeatMap(){};
	
		void getPic(OrderVector*);
		void getPic(SharedOrderVector*);

	private:
		void getSharedVectors(SharedOrderVector*);
		void sortSharedVectors();
		
		GlobalData* globaldata;
		vector<SharedRAbundVector*> lookup;
		RAbundVector rabund;
		string format, sorted, groupComb;
		ofstream outsvg;
		map<int, string> colorScale;
		map<int, string>::iterator it;

			
};
/***********************************************************************/

#endif



