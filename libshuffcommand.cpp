/*
 *  libshuffcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "libshuffcommand.h"

//**********************************************************************************************************************


LibShuffCommand::LibShuffCommand(){
	try {
		globaldata = GlobalData::getInstance();
		convert(globaldata->getCutOff(), cutOff);
		convert(globaldata->getIters(), iters);
		convert(globaldata->getStep(), step);
		form = globaldata->getForm();
		matrix = globaldata->gMatrix;
		coverageFile = getRootName(globaldata->getPhylipFile()) + "coverage";
		summaryFile = getRootName(globaldata->getPhylipFile()) + "slsummary";
		openOutputFile(coverageFile, out);
		openOutputFile(summaryFile, outSum);
		
		//set the groups to be analyzed
		setGroups();

		//file headers for coverage file
		out << "D" << '\t';
		for (int i = 0; i < groupComb.size(); i++) {
			out << "C" + groupComb[i] << '\t';
		}
		
		for (int i = 0; i < numGroups; i++) {
			for (int j = 0; j < numGroups; j++) {
				//don't output AA to AA
				if (i != j) {
					out << "Delta" + globaldata->Groups[i] + "-" + globaldata->Groups[j] << '\t';
				}
			}
		}
		out << endl;

		numComp = numGroups*numGroups;
		
		coverage = new Coverage();
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the LibShuffCommand class Function LibShuffCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the LibShuffCommand class function LibShuffCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
			
}

//**********************************************************************************************************************

LibShuffCommand::~LibShuffCommand(){
	delete coverage;
}

//**********************************************************************************************************************

int LibShuffCommand::execute(){
	try {
		//deltaValues[0] = scores for the difference between AA and AB.
		//cValues[0][0] = AA, cValues[0][1] = AB, cValues[0][2] = AC, cValues[1][0] = BA, cValues[1][1] = BB...
		vector<float> dist;
		int next;
		
		sumDelta.resize(numComp-numGroups, 0.0);
		
		float D = 0.0;
	
		/*****************************/
		//get values for users matrix
		/*****************************/
		matrix->setBounds();
		
		if (form != "discrete") { matrix->getDist(dist); next = 1; }
//cout << "Distances" << endl;
//for (int i = 0; i < dist.size(); i++) { cout << dist[i] << " "; }	
//cout << endl;
	
		//get values for users matrix
		while (D <= cutOff) {
			//clear out old Values
			deltaValues.clear();			
			coverage->getValues(matrix, D, cValues);
			
			//find delta values
			int count = 0;
			for (int i = 0; i < numGroups; i++) {
				for (int j = 0; j < numGroups; j++) {
					//don't save AA to AA
					if (i != j) {
						//(Caa - Cab)^2
						deltaValues.push_back( (abs(cValues[i][i]-cValues[i][j]) * abs(cValues[i][i]-cValues[i][j])) ); 
						sumDelta[count] += deltaValues[count];
						count++;
					}
				}
			}
			
			printCoverageFile(D);
			
			//check form
			if (form != "discrete") {   
				if (next == dist.size()) { break; }
				else {  D = dist[next];  next++;	}
			}else {  D += step;  }
			

		}
		
		//output sum Deltas
		for (int i = 0; i < numGroups; i++) {
			for (int j = 0; j < numGroups; j++) {
				//don't output AA to AA
				if (i != j) {
					cout << "Delta " + globaldata->Groups[i] + "-" + globaldata->Groups[j] << '\t';
				}
			}
		}
		cout << endl;
		
		for (int i = 0; i < sumDelta.size(); i++) {
			cout << setprecision(6) << sumDelta[i] << '\t';
		}
		cout << endl;
				
		/*******************************************************************************/
		//create and score random matrixes finding the sumDelta values for summary file
		/******************************************************************************/

		//initialize rsumDelta
		rsumDelta.resize(numComp-numGroups);
		for (int l = 0; l < rsumDelta.size(); l++) {
			for (int w = 0; w < iters; w++) {
				rsumDelta[l].push_back(0.0);
			}
		}
		
		
		for (int m = 0; m < iters; m++) {
			//generate random matrix in getValues
			//values for random matrix
			cout << "Iteration " << m+1 << endl;
			D = 0.0;
			next = 1;
			
			while (D <= cutOff) {
				coverage->getValues(matrix, D, cValues, "random");
			
				//find delta values
				int count = 0;
				for (int i = 0; i < numGroups; i++) {
					for (int j = 0; j < numGroups; j++) {
						//don't save AA to AA
						if (i != j) {
							//(Caa - Cab)^2
							rsumDelta[count][m] += ((abs(cValues[i][i]-cValues[i][j]) * abs(cValues[i][i]-cValues[i][j])));
//cout << "iter " << m << " box " << i << j << " delta = " << ((abs(cValues[i][i]-cValues[i][j]) * abs(cValues[i][i]-cValues[i][j]))) << endl;
							count++;
						}
					}
				}

				//check form
				if (form != "discrete") {   
					if (next == dist.size()) { break; }
					else {  D = dist[next];  next++;	}
				}else {  D += step;  }

			
				//clear out old Values
				cValues.clear();
			}
cout << "random sum delta for iter " << m << endl;
for (int i = 0; i < rsumDelta.size(); i++) {
	cout << setprecision(6) << rsumDelta[i][m] << '\t';
}
cout << endl;

		}
		
		/**********************************************************/
		//find the signifigance of the user matrix' sumdelta values
		/**********************************************************/
		
		for (int t = 0; t < rsumDelta.size(); t++) {
			//sort rsumDelta t
			sort(rsumDelta[t].begin(), rsumDelta[t].end());
			
			//the index of the score higher than yours is returned 
			//so if you have 1000 random matrices the index returned is 100 
			//then there are 900 matrices with a score greater then you. 
			//giving you a signifigance of 0.900
			int index = findIndex(sumDelta[t], t);    
			
			//the signifigance is the number of trees with the users score or higher 
			sumDeltaSig.push_back((iters-index)/(float)iters);

		}
		
		printSummaryFile();
		
		//clear out users groups
		globaldata->Groups.clear();
		
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the LibShuffCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the LibShuffCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************
void LibShuffCommand::printCoverageFile(float d) {
	try {
		//format output
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		out << setprecision(6) << d << '\t';
		
		//print out coverage values
		for (int i = 0; i < numGroups; i++) {
			for (int j = 0; j < numGroups; j++) {
				out << cValues[i][j] << '\t';
			}
		}
		
		//print out delta values
		for (int i = 0; i < deltaValues.size(); i++) {
			out << deltaValues[i] << '\t';
		}
		
		out << endl;
	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the LibShuffCommand class Function printCoverageFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the LibShuffCommand class function printCoverageFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
} 
//**********************************************************************************************************************
void LibShuffCommand::printSummaryFile() {
	try {
		//format output
		outSum.setf(ios::fixed, ios::floatfield); outSum.setf(ios::showpoint);
		
		for (int i = 0; i < numGroups; i++) {
			for (int j = 0; j < numGroups; j++) {
				//don't output AA to AA
				if (i != j) {
					outSum << "Delta " + globaldata->Groups[i] + "-" + globaldata->Groups[j] << '\t'<< "DeltaSig " + globaldata->Groups[i] + "-" + globaldata->Groups[j] << '\t';
					cout << "Delta " + globaldata->Groups[i] + "-" + globaldata->Groups[j] << '\t'<< "DeltaSig " + globaldata->Groups[i] + "-" + globaldata->Groups[j] << '\t';
				}
			}
		}
		outSum << endl;
		cout << endl;
		
		//print out delta values
		for (int i = 0; i < sumDelta.size(); i++) {
			outSum << setprecision(6) << sumDelta[i] << '\t' << setprecision(globaldata->getIters().length()) << sumDeltaSig[i] << '\t';
			cout << setprecision(6) << sumDelta[i] << '\t' << setprecision(globaldata->getIters().length()) << sumDeltaSig[i] << '\t';
		}
		
		outSum << endl;
		cout << endl;
	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the LibShuffCommand class Function printSummaryFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the LibShuffCommand class function printSummaryFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
} 

//**********************************************************************************************************************
void LibShuffCommand::setGroups() {
	try {
		//if the user has not entered specific groups to analyze then do them all
		if (globaldata->Groups.size() == 0) {
			numGroups = globaldata->gGroupmap->getNumGroups();
			for (int i=0; i < numGroups; i++) { 
				globaldata->Groups.push_back(globaldata->gGroupmap->namesOfGroups[i]);
			}
		}else {
			if (globaldata->getGroups() != "all") {
				//check that groups are valid
				for (int i = 0; i < globaldata->Groups.size(); i++) {
					if (globaldata->gGroupmap->isValidGroup(globaldata->Groups[i]) != true) {
						cout << globaldata->Groups[i] << " is not a valid group, and will be disregarded." << endl;
						// erase the invalid group from globaldata->Groups
						globaldata->Groups.erase (globaldata->Groups.begin()+i);
					}
				}
			
				//if the user only entered invalid groups
				if ((globaldata->Groups.size() == 0) || (globaldata->Groups.size() == 1)) { 
					numGroups = globaldata->gGroupmap->getNumGroups();
					for (int i=0; i < numGroups; i++) { 
						globaldata->Groups.push_back(globaldata->gGroupmap->namesOfGroups[i]);
					}
					cout << "When using the groups parameter you must have at least 2 valid groups. I will run the command using all the groups in your groupfile." << endl; 
				}else { numGroups = globaldata->Groups.size(); }
			}else { //users wants all groups
				numGroups = globaldata->gGroupmap->getNumGroups();
				globaldata->Groups.clear();
				for (int i=0; i < numGroups; i++) { 
					globaldata->Groups.push_back(globaldata->gGroupmap->namesOfGroups[i]);
				}
			}
		}
		
		//sort so labels match
		sort(globaldata->Groups.begin(), globaldata->Groups.end());
		
		// number of comparisons i.e. with groups A,B,C = AA, AB, AC, BA, BB, BC...;
		for (int i=0; i<numGroups; i++) { 
			for (int l = 0; l < numGroups; l++) {
				//set group comparison labels
				groupComb.push_back(globaldata->Groups[i] + "-" + globaldata->Groups[l]);
			}
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the LibShuffCommand class Function setGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the LibShuffCommand class function setGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/***********************************************************/
int LibShuffCommand::findIndex(float score, int index) {
	try{
		for (int i = 0; i < rsumDelta[index].size(); i++) {
			if (rsumDelta[index][i] >= score)	{	return i;	}
		}
		return rsumDelta[index].size();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the LibShuffCommand class Function findIndex. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the LibShuffCommand class function findIndex. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************/

