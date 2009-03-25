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
#include "sabundvector.hpp"
#include "sharedsabundvector.h"
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
		void setGroups();
	
		GlobalData* globaldata;
		vector<SharedSAbundVector> lookup;
		SAbundVector sabund;
		string format;

			
};

#endif

