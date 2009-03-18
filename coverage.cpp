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
void Coverage::getValues(FullMatrix* matrix, float d, vector< vector<float> >& data) {
	try {
		vector<float> min;
		vector<string> groups;
		
		//initialize data
		data.resize(numUserGroups);
		for (int l = 0; l < data.size(); l++) {
			data[l].push_back(0.0);
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
					
					min = matrix->getMins(count); //returns vector of mins for "box" requested ie. groups A, B, 0 = AA, 1 = AB, 2 = BA, 3 = BB;
					
					//find the coverage at this distance
					sort(min.begin(), min.end());
//cout << "minvector for : " << globaldata->gGroupmap->namesOfGroups[i] + globaldata->gGroupmap->namesOfGroups[j] << endl;
//for(int h = 0; h<min.size(); h++) {
// cout << min[h] << " ";
//}
//cout << endl;
					int index = -1;
					//find index in min where value is higher than d
					for (int m = 0; m < min.size(); m++) {
						if (min[m] > d)	{ index = m; break;	}
					}
					
					// if you don't find one than all the mins are less than d
					if (index == -1) { index = min.size();  }
					
					//save value in data
					data[a][b] = 1.0 - ((min.size()-index)/(float)min.size());
//cout << "D = " << d << "data " << a << b << " = " << data[a][b] << endl;		
					if (b < numUserGroups-1) {  b++;  }
					else{ //you are moving to a new row of "boxes"
						b = 0;
						a++;
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
//For the random matrices
void Coverage::getValues(FullMatrix* matrix, float d, vector< vector<float> >& data, string r) {
	try {
		vector<float> min;
		vector<string> groups;
		
		//initialize data
		data.resize(numUserGroups);
		for (int l = 0; l < data.size(); l++) {
			data[l].push_back(0.0);
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
cout << "combo " << a << b << endl;
cout << "original matrix mins4rows " << endl;
matrix->printMinsForRows(cout);					
					//create random matrix for this comparison
					matrix->shuffle(count);
			
					min = matrix->getMins(count);  //returns vector of mins for "box" requested ie. groups A, B, 0 = AA, 1 = AB, 2 = BA, 3 = BB;
cout << "shuffled matrix mins4rows " << endl;
matrix->printMinsForRows(cout);

					//find the coverage at this distance
					sort(min.begin(), min.end());
					
					int index = -1;
					//find index in min where value is higher than d
					for (int m = 0; m < min.size(); m++) {
						if (min[m] > d)	{ index = m; break;	}
					}
					
					// if you don't find one than all the mins are less than d
					if (index == -1) { index = min.size();  }
					
					//save value in data
					data[a][b] = 1.0 - ((min.size()-index)/(float)min.size());
cout << "D = " << d << "data " << a << b << " = " << data[a][b] << endl;		
					if (b < numUserGroups-1) {  b++;  }
					else{ //you are moving to a new row of "boxes"
						b = 0;
						a++;
					}
				}
				count++;
				
				//restore matrix to original form for next shuffle
				matrix->restore();
min = matrix->getMins(count-1); 
cout << "restored matrix mins4rows " << endl;
matrix->printMinsForRows(cout);
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

	