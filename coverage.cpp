/*
 *  coverage.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "coverage.h"

//**********************************************************************************************************************
vector< vector<float> > Coverage::getValues(FullMatrix*, float) {
	try {
		globaldata = GlobalData::getInstance();
		numGroups = globaldata->Groups.size();
		
		//initialize data
		data.resize(numGroups);
		for (int l = 0; l < data.size(); l++) {
			data[l].push_back(0.0);
		}

		/**************************************/
		//get the minumums for each comparision
		/**************************************/
		return data;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the Coverage class Function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the Coverage class function getValues. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************

	