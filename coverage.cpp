/*
 *  coverage.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class library coverage at the given distances of the Cramer-von Mises statistic.
	you may refer to the "Integration of Microbial Ecology and Statistics: A Test To Compare Gene Libraries" 
	paper in Applied and Environmental Microbiology, Sept. 2004, p. 5485-5492 0099-2240/04/$8.00+0  
	DOI: 10.1128/AEM.70.9.5485-5492.2004 Copyright 2004 American Society for Microbiology for more information. */
	
	
#include "coverage.h"

//**********************************************************************************************************************
Coverage::Coverage() {
		globaldata = GlobalData::getInstance();
		numUserGroups = globaldata->Groups.size();
		numGroups = globaldata->gGroupmap->getNumGroups();
}

//**********************************************************************************************************************
void Coverage::getValues(FullMatrix* matrix, vector< vector< vector<float> > >& data, vector<float> dist) {
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
		
		for (int i = 0; i < numGroups; i++) {
			for (int j = 0; j < numGroups; j++) {
			
				//is this "box" one the user wants analyzed?
				if ((inUsersGroups(globaldata->gGroupmap->namesOfGroups[i], globaldata->Groups) == true) && (inUsersGroups(globaldata->gGroupmap->namesOfGroups[j], globaldata->Groups) == true)) {
					
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
						data[k][i][j] = 1.0 - ((min.size()-index)/(float)min.size());
	
					}
				}
				count++;
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
void Coverage::getValues(FullMatrix* matrix, vector< vector< vector<float> > >& data, vector<float> dist, string mode) {
	try {
		vector<float> min1;
		vector<float> min2;
		vector<string> groups;
		
		//initialize data
		data.resize(dist.size());
		for (int l = 0; l < data.size(); l++) {
			data[l].resize(numGroups);
			for (int k = 0; k < data[l].size(); k++) {
				data[l][k].push_back(0.0);
			}
		}
		
		int count = 0;
		int count2 = 0;
		
		//for each box
		for (int i = 0; i < numGroups; i++) {
			for (int j = 0; j < numGroups; j++) {
				
				if (i != j) {
					//is this "box" one the user wants analyzed?
					if ((inUsersGroups(globaldata->gGroupmap->namesOfGroups[i], globaldata->Groups) == true) && (inUsersGroups(globaldata->gGroupmap->namesOfGroups[j], globaldata->Groups) == true)) {
					
						matrix->shuffle(globaldata->gGroupmap->namesOfGroups[i], globaldata->gGroupmap->namesOfGroups[j]);
					
						min1 = matrix->getMins(count);  //returns vector of mins for "box" requested ie. groups A, B, 0 = AA, 1 = AB, 2 = BA, 3 = BB;
						min2 = matrix->getMins(count2);  //returns vector of mins for "box" requested ie. groups A, B, 0 = AA, 1 = AB, 2 = BA, 3 = BB;

						//find the coverage at this distance
						sort(min1.begin(), min1.end());
					
						//find the coverage at this distance
						sort(min2.begin(), min2.end());
					
						float distDiff = 0;
						
						//loop through each distance and fill data
						for (int k = 0; k < data.size(); k++) {
							//****** coverage of AA **********//
							int index = -1;
							//find index in min where value is higher than d
							for (int m = 0; m < min1.size(); m++) {
									if (min1[m] > dist[k])	{ index = m; break;	}
							}
					
							// if you don't find one than all the mins are less than d
							if (index == -1) { index = min1.size();  }
							
							//****** coverage of AB **********//
							int index2 = -1;
							//find index in min where value is higher than d
							for (int m = 0; m < min2.size(); m++) {
									if (min2[m] > dist[k])	{ index2 = m; break;	}
							}
					
							// if you don't find one than all the mins are less than d
							if (index2 == -1) { index2 = min2.size();  }

							//coverage of ii
							float covII = 1.0 - ((min1.size()-index)/(float)min1.size());
							
							//coverage of ij
							float covIJ = 1.0 - ((min2.size()-index2)/(float)min2.size());
							
							//save value in data (Caa - Cab)^2 * distDiff
							data[k][i][j] = ((covII-covIJ) * (covII-covIJ)) * distDiff;
							
							//update distDiff
							if (k < data.size() - 1) {
								distDiff = dist[k+1] - dist[k];	
							}
						}
					
						//put matrix back to original
						matrix->restore();
					}
				}
				count2++;
			}
			count += numGroups+1; //go from AA to BB to CC
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

