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
Coverage::Coverage() {
		globaldata = GlobalData::getInstance();
		numUserGroups = globaldata->Groups.size();
		numGroups = globaldata->gGroupmap->getNumGroups();
}

//**********************************************************************************************************************
void Coverage::getValues(FullMatrix* matrix, vector< vector< vector<float> > >& data, vector<float> dist, string mode) {
	try {
		vector<float> min;
		vector<string> groups;
		
		//initialize data
		data.resize(dist.size());
		for (int l = 0; l < data.size(); l++) {
			data[l].resize(numGroups);
			for (int k = 0; k < data[l].size(); k++) {
				data[l][k].push_back(0.0);
			}
		}

		/**************************************/
		//get the minimums for each comparision
		/**************************************/
		int count = 0;
		int a = 0;
		int b = 0;
		for (int i = 0; i < numGroups; i++) {
			for (int j = 0; j < numGroups; j++) {
			
				//is this "box" one hte user wants analyzed?
				if ((inUsersGroups(globaldata->gGroupmap->namesOfGroups[i], globaldata->Groups) == true) && (inUsersGroups(globaldata->gGroupmap->namesOfGroups[j], globaldata->Groups) == true)) {
					
					if (mode == "random") {
						//create random matrix for this comparison
						matrix->shuffle(globaldata->gGroupmap->namesOfGroups[i], globaldata->gGroupmap->namesOfGroups[j]);
					}
			
					min = matrix->getMins(count);  //returns vector of mins for "box" requested ie. groups A, B, 0 = AA, 1 = AB, 2 = BA, 3 = BB;

					//find the coverage at this distance
					sort(min.begin(), min.end());
					
					//loop through each distance and fill data
					for (int k = 0; k < data.size(); k++) {
					
						int index = -1;
						//find index in min where value is higher than d
						for (int m = 0; m < min.size(); m++) {
							if (min[m] > dist[k])	{ index = m; break;	}
						}
					
						// if you don't find one than all the mins are less than d
						if (index == -1) { index = min.size();  }
					
						//save value in data
						data[k][a][b] = 1.0 - ((min.size()-index)/(float)min.size());
	
					}
					
					//move to next box
					if (b < numUserGroups-1) {  b++;  }
					else{ //you are moving to a new row of "boxes"
						b = 0;
						a++;
					}

					count++;
					
					if (mode == "random") {
						//restore matrix to original form for next shuffle
						matrix->restore();
					}
				}
			}
		}
		
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

	