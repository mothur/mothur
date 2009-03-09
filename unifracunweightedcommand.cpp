/*
 *  unifracunweightedcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 2/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "unifracunweightedcommand.h"

/***********************************************************/
UnifracUnweightedCommand::UnifracUnweightedCommand() {
	try {
		globaldata = GlobalData::getInstance();
		
		T = globaldata->gTree;
		tmap = globaldata->gTreemap;
		sumFile = globaldata->getTreeFile() + ".uwsummary";
		openOutputFile(sumFile, outSum);

		setGroups(); //sets users groups to analyze
		convert(globaldata->getIters(), iters);  //how many random trees to generate
		unweighted = new Unweighted(tmap);

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UnifracUnweightedCommand class Function UnifracUnweightedCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UnifracUnweightedCommand class function UnifracUnweightedCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/***********************************************************/
int UnifracUnweightedCommand::execute() {
	try {
	
		userData.resize(numComp,0);  //data[0] = unweightedscore 
		randomData.resize(numComp,0); //data[0] = unweightedscore
		//create new tree with same num nodes and leaves as users
				
		//get pscores for users trees
		for (int i = 0; i < T.size(); i++) {
			counter = 0;
			unweightedFile = globaldata->getTreeFile()  + toString(i+1) + ".unweighted";
			unweightedFileout = globaldata->getTreeFile() + "temp." + toString(i+1) + ".unweighted";
			
			//column headers
//			outSum << "Tree# " << i+1 << endl;
			outSum << "Tree#" << '\t' << "Groups" << '\t'  <<  "UWScore" <<'\t' << "UWSig" <<  endl;
//			cout << "Tree# " << i+1 << endl;
			cout << "Tree#" << '\t' << "Groups" << '\t'  <<  "UWScore" << '\t' << "UWSig" <<  endl;


			//get unweighted for users tree
			rscoreFreq.resize(numComp);  
			rCumul.resize(numComp);  
			utreeScores.resize(numComp);  
			UWScoreSig.resize(numComp); 

			userData = unweighted->getValues(T[i]);  //userData[0] = unweightedscore
			
			//output scores for each combination
			for(int k = 0; k < numComp; k++) {
				//saves users score
				utreeScores[k].push_back(userData[k]);
			}
			
			//get unweighted scores for random trees
			for (int j = 0; j < iters; j++) {
				//we need a different getValues because when we swap the labels we only want to swap those in each parwise comparison
				randomData = unweighted->getValues(T[i], "", "");
				
				for(int k = 0; k < numComp; k++) {	
					//add trees unweighted score to map of scores
					it2 = rscoreFreq[k].find(randomData[k]);
					if (it2 != rscoreFreq[k].end()) {//already have that score
						rscoreFreq[k][randomData[k]]++;
					}else{//first time we have seen this score
						rscoreFreq[k][randomData[k]] = 1;
					}
				
					//add randoms score to validscores
					validScores[randomData[k]] = randomData[k];
				}
			}
		
		for(int a = 0; a < numComp; a++) {
			float rcumul = 1.0000;
			//this loop fills the cumulative maps and put 0.0000 in the score freq map to make it easier to print.
			for (it = validScores.begin(); it != validScores.end(); it++) { 
				//make rscoreFreq map and rCumul
				it2 = rscoreFreq[a].find(it->first);
				rCumul[a][it->first] = rcumul;
				//get percentage of random trees with that info
				if (it2 != rscoreFreq[a].end()) {  rscoreFreq[a][it->first] /= iters; rcumul-= it2->second;  }
				else { rscoreFreq[a][it->first] = 0.0000; } //no random trees with that score
			}
			UWScoreSig[a].push_back(rCumul[a][userData[a]]);
		}
		
		printUnweightedFile();
		printUWSummaryFile();
		
		rscoreFreq.clear(); 
		rCumul.clear();  
		validScores.clear(); 
		utreeScores.clear();  
		UWScoreSig.clear(); 
	}
		//reset groups parameter
		globaldata->Groups.clear(); globaldata->setGroups("");
		outSum.close();
		
		return 0;
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UnifracUnweightedCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UnifracUnweightedCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/***********************************************************/
void UnifracUnweightedCommand::printUnweightedFile() {
	try {
		vector<double> data;
		
		for(int a = 0; a < numComp; a++) {
			initFile(groupComb[a]);
			//print each line
			for (it = validScores.begin(); it != validScores.end(); it++) { 
				data.push_back(it->first);  data.push_back(rscoreFreq[a][it->first]); data.push_back(rCumul[a][it->first]); 
				output(data);
				data.clear();
			} 
			resetFile();
		}
		
		out.close();
		inFile.close();
		remove(unweightedFileout.c_str());
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UnifracUnweightedCommand class Function printUnweightedFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UnifracUnweightedCommand class function printUnweightedFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************/
void UnifracUnweightedCommand::printUWSummaryFile() {
	try {
				
		//format output
		outSum.setf(ios::fixed, ios::floatfield); outSum.setf(ios::showpoint);
			
		//print each line

		for(int a = 0; a < numComp; a++) {
			if (UWScoreSig[a][0] > (1/(float)iters)) {
				outSum << setprecision(6) << groupComb[a] << '\t' << '\t' << utreeScores[a][0] << '\t' << setprecision(globaldata->getIters().length()) << UWScoreSig[a][0] << endl;
				cout << setprecision(6)  << groupComb[a] << '\t' << '\t' << utreeScores[a][0] << '\t' << setprecision(globaldata->getIters().length()) << UWScoreSig[a][0] << endl; 
			}else {
				outSum << setprecision(6) << groupComb[a] << '\t' << '\t' << utreeScores[a][0] << '\t' << setprecision(globaldata->getIters().length()) << "<" << (1/float(iters)) << endl;
				cout << setprecision(6)  << groupComb[a] << '\t' << '\t' << utreeScores[a][0] << '\t' << setprecision(globaldata->getIters().length()) << "<" << (1/float(iters)) << endl; 
			}
		}
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UnifracUnweightedCommand class Function printUWSummaryFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UnifracUnweightedCommand class function printUWSummaryFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/***********************************************************/

void UnifracUnweightedCommand::setGroups() {
	try {
		string allGroups = "";
		numGroups = 0;
		//if the user has not entered specific groups to analyze then do them all
		if (globaldata->Groups.size() != 0) {
			if (globaldata->Groups[0] != "all") {
				//check that groups are valid
				for (int i = 0; i < globaldata->Groups.size(); i++) {
					if (tmap->isValidGroup(globaldata->Groups[i]) != true) {
						cout << globaldata->Groups[i] << " is not a valid group, and will be disregarded." << endl;
						// erase the invalid group from globaldata->Groups
						globaldata->Groups.erase(globaldata->Groups.begin()+i);
					}
				}
			
				//if the user only entered invalid groups
				if (globaldata->Groups.size() == 0) { 
					cout << "When using the groups parameter you must have at least 1 valid group. I will run the command using all the groups in your groupfile." << endl; 
					for (int i = 0; i < tmap->namesOfGroups.size(); i++) {
						globaldata->Groups.push_back(tmap->namesOfGroups[i]);
						numGroups++;
						allGroups += tmap->namesOfGroups[i] + "-";
					}
					allGroups = allGroups.substr(0, allGroups.length()-1);
				}else {
					for (int i = 0; i < globaldata->Groups.size(); i++) {
						allGroups += globaldata->Groups[i] + "-";
						numGroups++;
					}
					allGroups = allGroups.substr(0, allGroups.length()-1);
				}
			}else{//user has enter "all" and wants the default groups
				globaldata->Groups.clear();
				for (int i = 0; i < tmap->namesOfGroups.size(); i++) {
					globaldata->Groups.push_back(tmap->namesOfGroups[i]);
					numGroups++;
					allGroups += tmap->namesOfGroups[i] + "-";
				}
				allGroups = allGroups.substr(0, allGroups.length()-1);
				globaldata->setGroups("");
			}
		}else {
			for (int i = 0; i < tmap->namesOfGroups.size(); i++) {
				allGroups += tmap->namesOfGroups[i] + "-";
			}
			allGroups = allGroups.substr(0, allGroups.length()-1);
			numGroups = 1;
		}
		
		//calculate number of comparsions
		numComp = 0;
		for (int r=0; r<numGroups; r++) { 
			for (int l = r+1; l < numGroups; l++) {
				groupComb.push_back(globaldata->Groups[r]+ "-" + globaldata->Groups[l]);
				numComp++;
			}
		}
		
		//ABC
		if (numComp != 1) {
			groupComb.push_back(allGroups);
			numComp++;
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UnifracUnweightedCommand class Function setGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UnifracUnweightedCommand class function setGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		

}
/*****************************************************************/

void UnifracUnweightedCommand::initFile(string label){
	try {
		if(counter != 0){
			openOutputFile(unweightedFileout, out);
			openInputFile(unweightedFile, inFile);

			string inputBuffer;
			getline(inFile, inputBuffer);
		
			out	<<  inputBuffer << '\t' << label + "RandFreq" << '\t' << label + "RandCumul" << endl;		
		}else{
			openOutputFile(unweightedFileout, out);
			out	<< label + "Score" << '\t' << label + "RandFreq" << '\t' << label + "RandCumul" << endl;
		}

		out.setf(ios::fixed, ios::floatfield);
		out.setf(ios::showpoint);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UnifracUnweightedCommand class Function initFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UnifracUnweightedCommand class function initFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

void UnifracUnweightedCommand::output(vector<double> data){
	try {
		if(counter != 0){		
			string inputBuffer;
			getline(inFile, inputBuffer);
//			out	<<  inputBuffer << setprecision(6) << '\t' << data[0] << setprecision(globaldata->getIters().length()) << '\t' << data[1] << '\t' << data[2] << endl;
		
			out << inputBuffer << '\t' << setprecision(6) << data[0] << setprecision(globaldata->getIters().length())  << '\t' << data[1] << '\t' << data[2] << endl;
		}
		else{
			out << setprecision(6) << data[0] << setprecision(globaldata->getIters().length())  << '\t' << data[1] << '\t' << data[2] << endl;
		}

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UnifracUnweightedCommand class Function output. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UnifracUnweightedCommand class function output. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
};

/***********************************************************************/

void UnifracUnweightedCommand::resetFile(){
	try {
		if(counter != 0){
			out.close();
			inFile.close();
		}
		else{
			out.close();
		}
		counter = 1;
		
		remove(unweightedFile.c_str());
		rename(unweightedFileout.c_str(), unweightedFile.c_str());
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UnifracUnweightedCommand class Function resetFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UnifracUnweightedCommand class function resetFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

