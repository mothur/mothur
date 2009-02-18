/*
 *  unifracweightedcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 2/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "unifracweightedcommand.h"

/***********************************************************/
UnifracWeightedCommand::UnifracWeightedCommand() {
	try {
		globaldata = GlobalData::getInstance();
		
		T = globaldata->gTree;
		tmap = globaldata->gTreemap;
		weightedFile = globaldata->getTreeFile() + ".weighted";
		openOutputFile(weightedFile, out);
		sumFile = globaldata->getTreeFile() + ".wsummary";
		openOutputFile(sumFile, outSum);
		distFile = globaldata->getTreeFile() + ".wdistrib";
		openOutputFile(distFile, outDist);
		
		//if the user has not entered specific groups to analyze then do them all
		if (globaldata->Groups.size() == 0) {
			numGroups = tmap->getNumGroups();
		}else {
			//check that groups are valid
			for (int i = 0; i < globaldata->Groups.size(); i++) {
				if (tmap->isValidGroup(globaldata->Groups[i]) != true) {
					cout << globaldata->Groups[i] << " is not a valid group, and will be disregarded." << endl;
					// erase the invalid group from globaldata->Groups
					globaldata->Groups.erase (globaldata->Groups.begin()+i);
				}
			}
			
			//if the user only entered invalid groups
			if (globaldata->Groups.size() == 0) { 
				numGroups = tmap->getNumGroups();
				cout << "When using the groups parameter you must have at least 2 valid groups. I will run the command using all the groups in your groupfile." << endl; 
			}else if (globaldata->Groups.size() == 1) { 
				cout << "When using the groups parameter you must have at least 2 valid groups. I will run the command using all the groups in your groupfile." << endl;
				numGroups = tmap->getNumGroups();
				globaldata->Groups.clear();
			}else { numGroups = globaldata->Groups.size(); }
		}
		
		//calculate number of comparisons i.e. with groups A,B,C = AB, AC, BC = 3;
		numComp = 0;
		int n = 1;
		for (int i=1; i<numGroups; i++) { 
			numComp += i; 
			for (int l = n; l < numGroups; l++) {
				//set group comparison labels
				if (globaldata->Groups.size() != 0) {
					groupComb.push_back(globaldata->Groups[i-1]+globaldata->Groups[l]);
				}else {
					groupComb.push_back(tmap->namesOfGroups[i-1]+tmap->namesOfGroups[l]);
				}
			}
			n++;
		}
			
		convert(globaldata->getIters(), iters);  //how many random trees to generate
		weighted = new Weighted(tmap);

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UnifracWeightedCommand class Function UnifracWeightedCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UnifracWeightedCommand class function UnifracWeightedCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/***********************************************************/
int UnifracWeightedCommand::execute() {
	try {
		
		//get weighted for users tree
		userData.resize(numComp,0);  //data[0] = weightedscore AB, data[1] = weightedscore AC...
		randomData.resize(numComp,0); //data[0] = weightedscore AB, data[1] = weightedscore AC...
		uscoreFreq.resize(numComp);  
		validScores.resize(numComp);  
		totalrscoreFreq.resize(numComp); 
		uCumul.resize(numComp); 	
		
		//format output
		outDist.setf(ios::fixed, ios::floatfield); outDist.setf(ios::showpoint);
		outDist << "Tree#" << '\t' << "Iter" << '\t' << "Groups"<< '\t' << "WScore" << endl;

		
		//create new tree with same num nodes and leaves as users
		randT = new Tree();
		
		//get pscores for users trees
		for (int i = 0; i < T.size(); i++) {
			rscoreFreq.resize(numComp);  //data[0] = weightedscore AB, data[1] = weightedscore AC...
			rCumul.resize(numComp); //data[0] = weightedscore AB, data[1] = weightedscore AC...	
				
			cout << "Processing tree " << i+1 << endl;
			userData = weighted->getValues(T[i]);  //userData[0] = weightedscore
			
			//save users score
			for (int s=0; s<numComp; s++) {
				//update uscoreFreq
				it = uscoreFreq[s].find(userData[s]);
				if (it == uscoreFreq[s].end()) {//new score
					uscoreFreq[s][userData[s]] = 1;
				}else{ uscoreFreq[s][userData[s]]++; }
				
				//add user score to valid scores
				validScores[s][userData[s]] = userData[s];

				//save users tree score for summary file
				utreeScores.push_back(userData[s]);
			}
			
			//copy T[i]'s info.
			randT->getCopy(T[i]); 
			
			//get scores for random trees
			for (int j = 0; j < iters; j++) {
				//create a random tree with same topology as T[i], but different labels
				randT->assembleRandomUnifracTree();
				//get pscore of random tree
				randomData = weighted->getValues(randT);
				
				//save ramdoms score
				for (int p=0; p<numComp; p++) {
					//add trees weighted score random score freq
					it2 = rscoreFreq[p].find(randomData[p]);
					if (it2 != rscoreFreq[p].end()) {//already have that score
						rscoreFreq[p][randomData[p]]++;
					}else{//first time we have seen this score
						rscoreFreq[p][randomData[p]] = 1;
					}
					
					//add random score to valid scores
					validScores[p][randomData[p]] = randomData[p];
					
					//output info to uwdistrib file
					outDist << i+1 << '\t' << '\t'<< j+1 << '\t' << '\t' << groupComb[p] << '\t'<< randomData[p] << endl;
				}
			}
			
			saveRandomScores(); //save all random scores for weighted file
			
			//find the signifigance of the score for summary file
			for (int t = 0; t < numComp; t++) {
				float rcumul = 0.0000;
				for (it = validScores[t].begin(); it != validScores[t].end(); it++) { 
					//make rscoreFreq map and rCumul
					it2 = rscoreFreq[t].find(it->first);
					//get percentage of random trees with that info
					if (it2 != rscoreFreq[t].end()) {  rscoreFreq[t][it->first] /= iters; rcumul+= it2->second;  }
					else { rscoreFreq[t][it->first] = 0.0000; } //no random trees with that score
					rCumul[t][it->first] = rcumul;
				}
			}
			
			//save the signifigance of the users score for printing later
			for (int f = 0; f < numComp; f++) {
				WScoreSig.push_back(rCumul[f][userData[f]]);
			}
			
			
			//clear random data
			rscoreFreq.clear();
			rCumul.clear();
		}
		
		rCumul.resize(numComp);
		for (int b = 0; b < numComp; b++) {
			float ucumul = 0.0000;
			float rcumul = 0.0000;
			//this loop fills the cumulative maps and put 0.0000 in the score freq map to make it easier to print.
			for (it = validScores[b].end(); it == validScores[b].begin(); it--) { 
				it2 = uscoreFreq[b].find(it->first);
				//user data has that score 
				if (it2 != uscoreFreq[b].end()) { uscoreFreq[b][it->first] /= T.size(); ucumul+= it2->second;  }
				else { uscoreFreq[b][it->first] = 0.0000; } //no user trees with that score
				//make uCumul map
				uCumul[b][it->first] = ucumul;
			
				//make rscoreFreq map and rCumul
				it2 = totalrscoreFreq[b].find(it->first);
				//get percentage of random trees with that info
				if (it2 != totalrscoreFreq[b].end()) {  totalrscoreFreq[b][it->first] /= (iters * T.size()); rcumul+= it2->second;  }
				else { totalrscoreFreq[b][it->first] = 0.0000; } //no random trees with that score
				rCumul[b][it->first] = rcumul;
			}
		}
		
		printWeightedFile();
		printWSummaryFile();
		
		//reset randomTree parameter to 0
		globaldata->setRandomTree("0");
		//clear out users groups
		globaldata->Groups.clear();
		
		delete randT;
		
		return 0;
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UnifracWeightedCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UnifracWeightedCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/***********************************************************/
void UnifracWeightedCommand::printWeightedFile() {
	try {
		//column headers
		
		out << "Group" << '\t' << "Score" << '\t' << "UserFreq" << '\t' << "UserCumul" << '\t' << "RandFreq" << '\t' << "RandCumul" << endl;
				
		//format output
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		//for each group
		for (int e = 0; e < numComp; e++) {
			//print each line in that group
			for (it = validScores[e].begin(); it != validScores[e].end(); it++) { 
				out << setprecision(6) <<  groupComb[e] << '\t' << it->first << '\t' << '\t' << uscoreFreq[e][it->first] << '\t' << uCumul[e][it->first] << '\t' << totalrscoreFreq[e][it->first] << '\t' << rCumul[e][it->first] << endl; 
			} 
		}
		
		out.close();
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UnifracWeightedCommand class Function printWeightedFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UnifracWeightedCommand class function printWeightedFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}


/***********************************************************/
void UnifracWeightedCommand::printWSummaryFile() {
	try {
		//column headers
		outSum << "Tree#" << '\t' << "Groups" << '\t' << '\t' << "WScore" << '\t' << '\t' << "WSig" <<  endl;
		
		//format output
		outSum.setf(ios::fixed, ios::floatfield); outSum.setf(ios::showpoint);
		
		//print each line
		int count = 0;
		for (int i = 0; i < T.size(); i++) { 
			for (int j = 0; j < numComp; j++) {
				outSum << setprecision(6) << i+1 << '\t' << '\t' << groupComb[j] << '\t' << utreeScores[count] << '\t' << WScoreSig[count] << endl; 
				count++;
			}
		}
		outSum.close();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UnifracWeightedCommand class Function printWeightedFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UnifracWeightedCommand class function printWeightedFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************/
void UnifracWeightedCommand::saveRandomScores() {
	try {
		for (int e = 0; e < numComp; e++) {
			//update total map with new random scores
			for (it = rscoreFreq[e].begin(); it != rscoreFreq[e].end(); it++) { 
				//does this score already exist in the total map
				it2 = totalrscoreFreq[e].find(it->first);
				//if yes then add them
				if (it2 != totalrscoreFreq[e].end()) { 
					totalrscoreFreq[e][it->first] += rscoreFreq[e][it->first];
				}else{ //its a new score
					totalrscoreFreq[e][it->first] = rscoreFreq[e][it->first];
				}
			}
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UnifracWeightedCommand class Function saveRandomScores. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UnifracWeightedCommand class function saveRandomScores. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************/
