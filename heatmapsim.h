#ifndef HEATMAPSIM_H
#define HEATMAPSIM_H
/*
 *  heatmapsim.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 6/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */


#include "sharedrabundvector.h"
#include "datavector.hpp"
#include "calculator.h"

/***********************************************************************/

class HeatMapSim {
	
	public:
		HeatMapSim(string, string, int);
		~HeatMapSim(){};
	
		vector<string> getPic(vector<SharedRAbundVector*>, vector<Calculator*>);
		string getPic(vector< vector<double> >, vector<string>);

	private:
		void printLegend(int, float);

		string format, groupComb, outputDir, inputfile;
		int fontSize;
		ofstream outsvg;
		MothurOut* m;
			
};

/***********************************************************************/

#endif

