/*
 *  libshuffcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 3/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

/* This class is designed to implement an integral form of the Cramer-von Mises statistic.
	you may refer to the "Integration of Microbial Ecology and Statistics: A Test To Compare Gene Libraries" 
	paper in Applied and Environmental Microbiology, Sept. 2004, p. 5485-5492 0099-2240/04/$8.00+0  
	DOI: 10.1128/AEM.70.9.5485-5492.2004 Copyright 2004 American Society for Microbiology for more information. */


#include "libshuffcommand.h"
#include "libshuff.h"
#include "slibshuff.h"
#include "dlibshuff.h"

//**********************************************************************************************************************

LibShuffCommand::LibShuffCommand(string option){
	try {
		srand( (unsigned)time( NULL ) );
		
		globaldata = GlobalData::getInstance();
		abort = false;
		Groups.clear();
		
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"iters","groups","step","form","cutoff"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			parser = new OptionParser();
			parser->parse(option, parameters);  delete parser;
			
			ValidParameters* validParameter = new ValidParameters();
		
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter->isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//make sure the user has already run the read.dist command
			if ((globaldata->gMatrix == NULL) || (globaldata->gGroupmap == NULL)) {
				cout << "You must read in a matrix and groupfile using the read.dist command, before you use the libshuff command. " << endl; abort = true;; 
			}
						
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			groups = validParameter->validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; savegroups = groups; }
			else { 
				savegroups = groups;
				splitAtDash(groups, Groups);
				globaldata->Groups = Groups;
			}
				
			string temp;
			temp = validParameter->validFile(parameters, "iters", false);				if (temp == "not found") { temp = "10000"; }
			convert(temp, iters); 
			
			temp = validParameter->validFile(parameters, "cutoff", false);				if (temp == "not found") { temp = "1.0"; }
			convert(temp, cutOff); 
			
			temp = validParameter->validFile(parameters, "step", false);				if (temp == "not found") { temp = "0.01"; }
			convert(temp, step); 
	
			userform = validParameter->validFile(parameters, "form", false);			if (userform == "not found") { userform = "integral"; }
			
			delete validParameter;
			
			if (abort == false) {
			
				matrix = globaldata->gMatrix;				//get the distance matrix
				setGroups();								//set the groups to be analyzed

				if(userform == "discrete"){
					form = new DLibshuff(matrix, iters, step, cutOff);
				}
				else{
					form = new SLibshuff(matrix, iters, cutOff);
				}
			}
			
		}
		
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

void LibShuffCommand::help(){
	try {
		cout << "The libshuff command can only be executed after a successful read.dist command including a groupfile." << "\n";
		cout << "The libshuff command parameters are groups, iters, step, form and cutoff.  No parameters are required." << "\n";
		cout << "The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 2 valid groups." << "\n";
		cout << "The group names are separated by dashes.  The iters parameter allows you to specify how many random matrices you would like compared to your matrix." << "\n";
		cout << "The step parameter allows you to specify change in distance you would like between each output if you are using the discrete form." << "\n";
		cout << "The form parameter allows you to specify if you would like to analyze your matrix using the discrete or integral form. Your options are integral or discrete." << "\n";
		cout << "The libshuff command should be in the following format: libshuff(groups=yourGroups, iters=yourIters, cutOff=yourCutOff, form=yourForm, step=yourStep)." << "\n";
		cout << "Example libshuff(groups=A-B-C, iters=500, form=discrete, step=0.01, cutOff=2.0)." << "\n";
		cout << "The default value for groups is all the groups in your groupfile, iters is 10000, cutoff is 1.0, form is integral and step is 0.01." << "\n";
		cout << "The libshuff command output two files: .coverage and .slsummary their descriptions are in the manual." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. iters), '=' and parameters (i.e.yourIters)." << "\n" << "\n";
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the LibShuffCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the LibShuffCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************

int LibShuffCommand::execute(){
	try {
		
		if (abort == true) {	return 0;	}
		
		savedDXYValues = form->evaluateAll();
		savedMinValues = form->getSavedMins();
		
		pValueCounts.resize(numGroups);
		for(int i=0;i<numGroups;i++){
			pValueCounts[i].assign(numGroups, 0);
		}
		
		Progress* reading = new Progress();
		
		for(int i=0;i<numGroups-1;i++) {
			for(int j=i+1;j<numGroups;j++) {
				reading->newLine(groupNames[i]+'-'+groupNames[j], iters);
				for(int p=0;p<iters;p++) {		
					form->randomizeGroups(i,j);
					if(form->evaluatePair(i,j) >= savedDXYValues[i][j])	{	pValueCounts[i][j]++;	}
					if(form->evaluatePair(j,i) >= savedDXYValues[j][i])	{	pValueCounts[j][i]++;	}
					reading->update(p);			
				}
				form->resetGroup(i);
				form->resetGroup(j);
			}
		}
		reading->finish();
		delete reading;

		cout << endl;
		printSummaryFile();
		printCoverageFile();
		
		//clear out users groups
		globaldata->Groups.clear();
		delete form;
		
		//delete globaldata's copy of the gmatrix to free up memory
		delete globaldata->gMatrix;  globaldata->gMatrix = NULL;
		
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

		ofstream outCov;
		summaryFile = getRootName(globaldata->getPhylipFile()) + "libshuff.coverage";
		openOutputFile(summaryFile, outCov);
		outCov.setf(ios::fixed, ios::floatfield); outCov.setf(ios::showpoint);
		cout.setf(ios::fixed, ios::floatfield); cout.setf(ios::showpoint);
		
		map<double,vector<int> > allDistances;
		map<double,vector<int> >::iterator it;

		vector<vector<int> > indices(numGroups);
		int numIndices = numGroups * numGroups;
		
		int index = 0;
		for(int i=0;i<numGroups;i++){
			indices[i].assign(numGroups,0);
			for(int j=0;j<numGroups;j++){
				indices[i][j] = index++;
				for(int k=0;k<savedMinValues[i][j].size();k++){
					if(allDistances[savedMinValues[i][j][k]].size() != 0){
						allDistances[savedMinValues[i][j][k]][indices[i][j]]++;
					}
					else{
						allDistances[savedMinValues[i][j][k]].assign(numIndices, 0);
						allDistances[savedMinValues[i][j][k]][indices[i][j]] = 1;
					}
				}
			}
		}
		it=allDistances.begin();
		
		cout << setprecision(8);

		vector<int> prevRow = it->second;
		it++;
		
		for(it;it!=allDistances.end();it++){
			for(int i=0;i<it->second.size();i++){
				it->second[i] += prevRow[i];
			}
			prevRow = it->second;
		}
		
		vector<int> lastRow = allDistances.rbegin()->second;
		outCov << setprecision(8);
		
		outCov << "dist";
		for (int i = 0; i < numGroups; i++){
			outCov << '\t' << groupNames[i];
		}
		for (int i=0;i<numGroups;i++){
			for(int j=i+1;j<numGroups;j++){
				outCov << '\t' << groupNames[i] << '-' << groupNames[j] << '\t';
				outCov << groupNames[j] << '-' << groupNames[i];
			}
		}
		outCov << endl;
		
		for(it=allDistances.begin();it!=allDistances.end();it++){
			outCov << it->first << '\t';
			for(int i=0;i<numGroups;i++){
				outCov << it->second[indices[i][i]]/(float)lastRow[indices[i][i]] << '\t';
			}
			for(int i=0;i<numGroups;i++){
				for(int j=i+1;j<numGroups;j++){
					outCov << it->second[indices[i][j]]/(float)lastRow[indices[i][j]] << '\t';
					outCov << it->second[indices[j][i]]/(float)lastRow[indices[j][i]] << '\t';
				}
			}
			outCov << endl;
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

		ofstream outSum;
		summaryFile = getRootName(globaldata->getPhylipFile()) + "libshuff.summary";
		openOutputFile(summaryFile, outSum);

		outSum.setf(ios::fixed, ios::floatfield); outSum.setf(ios::showpoint);
		cout.setf(ios::fixed, ios::floatfield); cout.setf(ios::showpoint);
		
		cout << setw(20) << left << "Comparison" << '\t' << setprecision(8) << "dCXYScore" << '\t' << "Significance" << endl;
		outSum << setw(20) << left << "Comparison" << '\t' << setprecision(8) << "dCXYScore" << '\t' << "Significance" << endl;
	
		int precision = (int)log10(iters);
		for(int i=0;i<numGroups;i++){
			for(int j=i+1;j<numGroups;j++){
				if(pValueCounts[i][j]){
					cout << setw(20) << left << groupNames[i]+'-'+groupNames[j] << '\t' << setprecision(8) << savedDXYValues[i][j] << '\t' << setprecision(precision) << pValueCounts[i][j]/(float)iters << endl;
					outSum << setw(20) << left << groupNames[i]+'-'+groupNames[j] << '\t' << setprecision(8) << savedDXYValues[i][j] << '\t' << setprecision(precision) << pValueCounts[i][j]/(float)iters << endl;
				}
				else{
					cout << setw(20) << left << groupNames[i]+'-'+groupNames[j] << '\t' << setprecision(8) << savedDXYValues[i][j] << '\t' << '<' <<setprecision(precision) << 1/(float)iters << endl;
					outSum << setw(20) << left << groupNames[i]+'-'+groupNames[j] << '\t' << setprecision(8) << savedDXYValues[i][j] << '\t' << '<' <<setprecision(precision) << 1/(float)iters << endl;
				}
				if(pValueCounts[j][i]){
					cout << setw(20) << left << groupNames[j]+'-'+groupNames[i] << '\t' << setprecision(8) << savedDXYValues[j][i] << '\t' << setprecision (precision) << pValueCounts[j][i]/(float)iters << endl;
					outSum << setw(20) << left << groupNames[j]+'-'+groupNames[i] << '\t' << setprecision(8) << savedDXYValues[j][i] << '\t' << setprecision (precision) << pValueCounts[j][i]/(float)iters << endl;
				}
				else{
					cout << setw(20) << left << groupNames[j]+'-'+groupNames[i] << '\t' << setprecision(8) << savedDXYValues[j][i] << '\t' << '<' <<setprecision (precision) << 1/(float)iters << endl;
					outSum << setw(20) << left << groupNames[j]+'-'+groupNames[i] << '\t' << setprecision(8) << savedDXYValues[j][i] << '\t' << '<' <<setprecision (precision) << 1/(float)iters << endl;
				}
			}
		}
		
		
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
		} else {
			if (savegroups != "all") {
				//check that groups are valid
				for (int i = 0; i < globaldata->Groups.size(); i++) {
					if (globaldata->gGroupmap->isValidGroup(globaldata->Groups[i]) != true) {
						cout << globaldata->Groups[i] << " is not a valid group, and will be disregarded." << endl;
						// erase the invalid group from globaldata->Groups
						globaldata->Groups.erase(globaldata->Groups.begin()+i);
					}
				}
			
				//if the user only entered invalid groups
				if ((globaldata->Groups.size() == 0) || (globaldata->Groups.size() == 1)) { 
					numGroups = globaldata->gGroupmap->getNumGroups();
					for (int i=0; i < numGroups; i++) { 
						globaldata->Groups.push_back(globaldata->gGroupmap->namesOfGroups[i]);
					}
					cout << "When using the groups parameter you must have at least 2 valid groups. I will run the command using all the groups in your groupfile." << endl; 
				} else { numGroups = globaldata->Groups.size(); }
			} else { //users wants all groups
				numGroups = globaldata->gGroupmap->getNumGroups();
				globaldata->Groups.clear();
				for (int i=0; i < numGroups; i++) { 
					globaldata->Groups.push_back(globaldata->gGroupmap->namesOfGroups[i]);
				}
			}
		}

		//sort so labels match
		sort(globaldata->Groups.begin(), globaldata->Groups.end());
		
		//sort
		sort(globaldata->gGroupmap->namesOfGroups.begin(), globaldata->gGroupmap->namesOfGroups.end());

		groupNames = globaldata->Groups;

		// number of comparisons i.e. with groups A,B,C = AA, AB, AC, BA, BB, BC...;
//		for (int i=0; i<numGroups; i++) { 
//			for (int l = 0; l < numGroups; l++) {
//				//set group comparison labels
//				groupComb.push_back(globaldata->Groups[i] + "-" + globaldata->Groups[l]);
//			}
//		}
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
