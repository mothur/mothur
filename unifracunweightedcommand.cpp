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
		unweightedFile = globaldata->getTreeFile() + ".unweighted";
		openOutputFile(unweightedFile, out);
		//column headers
		out << "Comb" << '\t' << "Score" << '\t' << "UserFreq" << '\t' << "UserCumul" << '\t' << "RandFreq" << '\t' << "RandCumul" << endl;
				
		sumFile = globaldata->getTreeFile() + ".uwsummary";
		openOutputFile(sumFile, outSum);
		//column headers
		outSum << "Tree#" << '\t' << "Comb" << '\t'  <<  "UWScore" << '\t' << '\t' << "UWSig" <<  endl;

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
		randT = new Tree();
				
		//get pscores for users trees
		for (int i = 0; i < T.size(); i++) {
			//get unweighted for users tree
			rscoreFreq.resize(numComp);  
			uscoreFreq.resize(numComp);  
			rCumul.resize(numComp);  
			uCumul.resize(numComp);  
			validScores.resize(numComp); 
			utreeScores.resize(numComp);  
			UWScoreSig.resize(numComp); 

			cout << "Processing tree " << i+1 << endl;
			outSum << "Tree#" << i+1 << endl;
			out << "Tree#" << i+1 << endl;
			userData = unweighted->getValues(T[i]);  //userData[0] = unweightedscore
			
			//output scores for each combination
			for(int k = 0; k < numComp; k++) {
				//update uscoreFreq
				it = uscoreFreq[k].find(userData[k]);
				if (it == uscoreFreq[k].end()) {//new score
					uscoreFreq[k][userData[k]] = 1;
				}else{ uscoreFreq[k][userData[k]]++; }
			
				//add users score to valid scores
				validScores[k][userData[k]] = userData[k];
			
				//saves users score
				utreeScores[k].push_back(userData[k]);
			}
			
			//copy T[i]'s info.
			randT->getCopy(T[i]); 
			
			//get unweighted scores for random trees
			for (int j = 0; j < iters; j++) {
				//we need a different getValues because when we swap the labels we only want to swap those in each parwise comparison
				randomData = unweighted->getValues(randT, "", "");
				
				for(int k = 0; k < numComp; k++) {	
cout << "iter " << j << " comp " << k << " = " << randomData[k] << endl;
					//add trees unweighted score to map of scores
					it2 = rscoreFreq[k].find(randomData[k]);
					if (it2 != rscoreFreq[k].end()) {//already have that score
						rscoreFreq[k][randomData[k]]++;
					}else{//first time we have seen this score
						rscoreFreq[k][randomData[k]] = 1;
					}
				
					//add randoms score to validscores
					validScores[k][randomData[k]] = randomData[k];
				}
			}
		
		for(int a = 0; a < numComp; a++) {
			float ucumul = 1.0000;
			float rcumul = 1.0000;
			//this loop fills the cumulative maps and put 0.0000 in the score freq map to make it easier to print.
			for (it = validScores[a].begin(); it != validScores[a].end(); it++) { 
				it2 = uscoreFreq[a].find(it->first);
				//make uCumul map
				uCumul[a][it->first] = ucumul;
				//user data has that score 
				if (it2 != uscoreFreq[a].end()) { uscoreFreq[a][it->first] /= T.size(); ucumul-= it2->second;  }
				else { uscoreFreq[a][it->first] = 0.0000; } //no user trees with that score
						
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
		uscoreFreq.clear();  
		rCumul.clear();  
		uCumul.clear();  
		validScores.clear(); 
		utreeScores.clear();  
		UWScoreSig.clear(); 
	}
		//reset groups parameter
		globaldata->Groups.clear(); globaldata->setGroups("");
		
		delete randT;
		
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
		//format output
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		for(int a = 0; a < numComp; a++) {
			//print each line
			for (it = validScores[a].begin(); it != validScores[a].end(); it++) { 
				out << setprecision(6) << groupComb[a] << '\t' << it->first << '\t' << '\t' << uscoreFreq[a][it->first] << '\t' << uCumul[a][it->first] << '\t' << rscoreFreq[a][it->first] << '\t' << rCumul[a][it->first] << endl; 
			} 
		}
		out.close();
		
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
		for (int i = 0; i< T.size(); i++) {
			for(int a = 0; a < numComp; a++) {
				outSum << setprecision(6) << i+1 << '\t' << groupComb[a] << '\t' << '\t' << utreeScores[a][i] << '\t' << UWScoreSig[a][i] << endl;
				cout << setprecision(6) << i+1 << '\t' << groupComb[a] << '\t' << '\t' << utreeScores[a][i] << '\t' << UWScoreSig[a][i] << endl; 
			}	
		}
		
		outSum.close();
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
						allGroups += tmap->namesOfGroups[i];
					}
				}else {
					for (int i = 0; i < globaldata->Groups.size(); i++) {
						allGroups += globaldata->Groups[i];
						numGroups++;
					}
				}
			}else{//user has enter "all" and wants the default groups
				for (int i = 0; i < tmap->namesOfGroups.size(); i++) {
					globaldata->Groups.push_back(tmap->namesOfGroups[i]);
					numGroups++;
					allGroups += tmap->namesOfGroups[i];
				}
				globaldata->setGroups("");
			}
		}else {
			for (int i = 0; i < tmap->namesOfGroups.size(); i++) {
				allGroups += tmap->namesOfGroups[i];
			}
			numGroups = 1;
		}
		
		//calculate number of comparsions
		numComp = 0;
		for (int r=0; r<numGroups; r++) { 
			for (int l = r+1; l < numGroups; l++) {
				groupComb.push_back(globaldata->Groups[r]+globaldata->Groups[l]);
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

