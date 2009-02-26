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
		//column headers
		out << "Group" << '\t' << "Score" << '\t' << "UserFreq" << '\t' << "UserCumul" << '\t' << "RandFreq" << '\t' << "RandCumul" << endl;

		sumFile = globaldata->getTreeFile() + ".wsummary";
		openOutputFile(sumFile, outSum);
				
		setGroups();	//sets the groups the user wants to analyze			
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
				
		//create new tree with same num nodes and leaves as users
		randT = new Tree();
		
		//get pscores for users trees
		for (int i = 0; i < T.size(); i++) {
			rScores.resize(numComp);  //data[0] = weightedscore AB, data[1] = weightedscore AC...
			uScores.resize(numComp);  //data[0] = weightedscore AB, data[1] = weightedscore AC...
			validScores.resize(numComp); 
							
			cout << "Processing tree " << i+1 << endl;
			userData = weighted->getValues(T[i]);  //userData[0] = weightedscore
			
			//save users score
			for (int s=0; s<numComp; s++) {
				//add users score to vector of user scores
				uScores[s].push_back(userData[s]);
				
				//add users score to vector of valid scores
				validScores[s].push_back(userData[s]);

				//save users tree score for summary file
				utreeScores.push_back(userData[s]);
			}
			
			//get scores for random trees
			for (int j = 0; j < iters; j++) {
//				int n = 1;
				int count = 0;
				for (int r=0; r<numGroups; r++) { 
					for (int l = r+1; l < numGroups; l++) {
						//copy T[i]'s info.
						randT->getCopy(T[i]);
						 
						if (globaldata->Groups.size() != 0) {
							//create a random tree with same topology as T[i], but different labels
							randT->assembleRandomUnifracTree(globaldata->Groups[r], globaldata->Groups[l]);
							//get wscore of random tree
							randomData = weighted->getValues(randT, globaldata->Groups[r], globaldata->Groups[l]);
						}else {
							//create a random tree with same topology as T[i], but different labels
							randT->assembleRandomUnifracTree(tmap->namesOfGroups[r], tmap->namesOfGroups[l]);
							//get wscore of random tree
							randomData = weighted->getValues(randT, tmap->namesOfGroups[r], tmap->namesOfGroups[l]);
						}
//						randT->createNewickFile("hold"+toString(r)+toString(l)+toString(j));

						//save scores
						rScores[count].push_back(randomData[0]);
						validScores[count][randomData[0]] = randomData[0];
						count++;
					}
//					n++;
				}
			}

			removeValidScoresDuplicates(); 
			//find the signifigance of the score for summary file
			for (int f = 0; f < numComp; f++) {
				//sort random scores
				sort(rScores[f].begin(), rScores[f].end());
				
				//the index of the score higher than yours is returned 
				//so if you have 1000 random trees the index returned is 100 
				//then there are 900 trees with a score greater then you. 
				//giving you a signifigance of 0.900
				int index = findIndex(userData[f], f);    if (index == -1) { cout << "error in UnifracWeightedCommand" << endl; exit(1); } //error code
			
				//the signifigance is the number of trees with the users score or higher 
				WScoreSig.push_back((iters-index)/(float)iters);
			}
			
			out << "Tree# " << i << endl;
			//printWeightedFile();
			
			//clear data
			rScores.clear();
			uScores.clear();
			validScores.clear();
		}
		
		printWSummaryFile();
		
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
/***********************************************************
void UnifracWeightedCommand::printWeightedFile() {
	try {
						
		//format output
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		//for each group
		for (int e = 0; e < numComp; e++) {
			//print each line in that group
			for (i = 0; i < validScores[e].size(); i++) { 
				out << setprecision(6) <<  groupComb[e] << '\t' << validScores[e][i] << '\t' << '\t' << uscoreFreq[e][it->first] << '\t' << uCumul[e][it->first] << '\t' << rscoreFreq[e][it->first] << '\t' << rCumul[e][it->first] << endl; 
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
		cout << "Tree#" << '\t' << "Groups" << '\t' << '\t' << "WScore" << '\t' << '\t' << "WSig" <<  endl;
		
		//format output
		outSum.setf(ios::fixed, ios::floatfield); outSum.setf(ios::showpoint);
		
		//print each line
		int count = 0;
		for (int i = 0; i < T.size(); i++) { 
			for (int j = 0; j < numComp; j++) {
				outSum << setprecision(6) << i+1 << '\t' << '\t' << groupComb[j] << '\t' << utreeScores[count] << '\t' << WScoreSig[count] << endl; 
				cout << setprecision(6) << i+1 << '\t' << '\t' << groupComb[j] << '\t' << utreeScores[count] << '\t' << WScoreSig[count] << endl; 
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
void UnifracWeightedCommand::removeValidScoresDuplicates() {
	try {
		for (int e = 0; e < numComp; e++) {
			//sort valid scores
			sort(validScores[e].begin(), validScores[e].end());
			
			for (int i = 0; i< validScores[e].size()-1; i++) { 
				if (validScores[e][i] == validScores[e][i+1]) { validScores[e].erase(validScores[e].begin()+i); }
			}
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UnifracWeightedCommand class Function removeValidScoresDuplicates. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UnifracWeightedCommand class function removeValidScoresDuplicates. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************/
int UnifracWeightedCommand::findIndex(float score, int index) {
	try{
		for (int i = 0; i < rScores[index].size(); i++) {
			if (rScores[index][i] >= score)	{	return i;	}
		}
		return rScores[index].size();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UnifracWeightedCommand class Function findIndex. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UnifracWeightedCommand class function findIndex. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************/
void UnifracWeightedCommand::setGroups() {
	try {
		//if the user has not entered specific groups to analyze then do them all
		if (globaldata->Groups.size() == 0) {
			numGroups = tmap->getNumGroups();
		}else {
			if (globaldata->getGroups() != "all") {
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
			}else { //users wants all groups
				numGroups = tmap->getNumGroups();
				globaldata->Groups.clear();
				globaldata->setGroups("");
			}
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
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UnifracWeightedCommand class Function setGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UnifracWeightedCommand class function setGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

