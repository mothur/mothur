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
		//cValues[0][0][0] = AA at distance 0.0, cValues[0][0][1] = AB at distance 0.0, cValues[0][0][2] = AC at distance 0.0, cValues[0][1][0] = BA at distance 0.0, cValues[0][1][1] = BB...
		Progress* reading;
		reading = new Progress("Comparing to random:", iters);
		
		sumDelta.resize(numComp-numGroups, 0.0);
		
		matrix->setBounds();
		
		//load distances
		if (form != "discrete") { matrix->getDist(dist); }
		else {
			float f = 0.0;
			while (f <= cutOff) {
				dist.push_back(f);
				f += step;
			}
		}
	
		/*****************************/
		//get values for users matrix
		/*****************************/
			
		//clear out old Values
		deltaValues.clear();
		deltaValues.resize(dist.size());			
		
		coverage->getValues(matrix, cValues, dist, "user");
		
		//loop through each distance and load rsumdelta
		for (int p = 0; p < cValues.size(); p++) {	
			//find delta values
			int count = 0;
			for (int i = 0; i < numGroups; i++) {
				for (int j = 0; j < numGroups; j++) {
					//don't save AA to AA
					if (i != j) {
						//(Caa - Cab)^2
						deltaValues[p].push_back( (abs(cValues[p][i][i]-cValues[p][i][j]) * abs(cValues[p][i][i]-cValues[p][i][j])) ); 
						sumDelta[count] += deltaValues[p][count];
						count++;
					}
				}
			}
		}
			
		printCoverageFile();
			
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
		
			coverage->getValues(matrix, cValues, dist, "random");
			
			//loop through each distance and load rsumdelta
			for (int p = 0; p < cValues.size(); p++) {
				//find delta values
				int count = 0;
				for (int i = 0; i < numGroups; i++) {
					for (int j = 0; j < numGroups; j++) {
						//don't save AA to AA
						if (i != j) {
							//(Caa - Cab)^2
							rsumDelta[count][m] += ((abs(cValues[p][i][i]-cValues[p][i][j]) * abs(cValues[p][i][i]-cValues[p][i][j])));
							count++;
						}
					}
				}
				
			}
//cout << "iter " << m << endl;
			//clear out old Values
			reading->update(m);
			cValues.clear();
			
//cout << "random sum delta for iter " << m << endl;
//for (int i = 0; i < rsumDelta.size(); i++) {
//	cout << setprecision(6) << rsumDelta[i][m] << '\t';
//}
//cout << endl;

		}
		
		reading->finish();
		delete reading;
				
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
void LibShuffCommand::printCoverageFile() {
	try {
		//format output
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		//loop through each distance 
		for (int p = 0; p < cValues.size(); p++) {
			out << setprecision(6) << dist[p] << '\t';
			//print out coverage values
			for (int i = 0; i < numGroups; i++) {
				for (int j = 0; j < numGroups; j++) {
					out << cValues[p][i][j] << '\t';
				}
			}
			
			for (int h = 0; h < deltaValues[p].size(); h++) {
				out << deltaValues[p][h] << '\t';
			}
			
			out << endl;
		}
		
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
			if (sumDeltaSig[i] > (1/(float)iters)) {
				outSum << setprecision(6) << sumDelta[i] << '\t' << setprecision(globaldata->getIters().length()) << sumDeltaSig[i] << '\t';
				cout << setprecision(6) << sumDelta[i] << '\t' << setprecision(globaldata->getIters().length()) << sumDeltaSig[i] << '\t';
			}else {
				outSum << setprecision(6) << sumDelta[i] << '\t' << setprecision(globaldata->getIters().length()) << "<" << (1/float(iters)) << '\t';
				cout << setprecision(6) << sumDelta[i] << '\t' << setprecision(globaldata->getIters().length()) << "<" << (1/float(iters)) << '\t';
			}
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

