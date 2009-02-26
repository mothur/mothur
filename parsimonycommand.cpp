/*
 *  parsimonycommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/26/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "parsimonycommand.h"

/***********************************************************/
ParsimonyCommand::ParsimonyCommand() {
	try {
		globaldata = GlobalData::getInstance();
		
		//randomtree will tell us if user had their own treefile or if they just want the random distribution
		randomtree = globaldata->getRandomTree();
		
		//user has entered their own tree
		if (randomtree == "") { 
			T = globaldata->gTree;
			tmap = globaldata->gTreemap;
			parsFile = globaldata->getTreeFile() + ".parsimony";
			openOutputFile(parsFile, out);
			sumFile = globaldata->getTreeFile() + ".psummary";
			openOutputFile(sumFile, outSum);
			//set users groups to analyze
			setGroups();
		}else { //user wants random distribution
			savetmap = globaldata->gTreemap;
			getUserInput();
			parsFile = randomtree + ".rd_parsimony";
			openOutputFile(parsFile, out);
		}
		
		convert(globaldata->getIters(), iters);  //how many random trees to generate
		pars = new Parsimony(tmap);

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ParsimonyCommand class Function ParsimonyCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ParsimonyCommand class function ParsimonyCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/***********************************************************/
int ParsimonyCommand::execute() {
	try {
		
				
		if (randomtree == "") {
			//get pscore for users tree
			userData.resize(numComp,0);  //data = AB, AC, BC, ABC.
			randomData.resize(numComp,0);  //data = AB, AC, BC, ABC.
			rscoreFreq.resize(numComp);  
			uscoreFreq.resize(numComp);  
			rCumul.resize(numComp);  
			uCumul.resize(numComp);  
			validScores.resize(numComp); 
			userTreeScores.resize(numComp);  
			UScoreSig.resize(numComp); 

			//get pscores for users trees
			for (int i = 0; i < T.size(); i++) {
				cout << "Processing tree " << i+1 << endl;
				userData = pars->getValues(T[i]);  //data = AB, AC, BC, ABC.
				
				//output scores for each combination
				for(int k = 0; k < numComp; k++) {
					cout << "Tree " << i+1 << " Combination " << groupComb[k] << " parsimony score = " << userData[k] << endl;
					//update uscoreFreq
					it = uscoreFreq[k].find(userData[k]);
					if (it == uscoreFreq[k].end()) {//new score
						uscoreFreq[k][userData[k]] = 1;
					}else{ uscoreFreq[k][userData[k]]++; }
					
					//add users score to valid scores
					validScores[k][userData[k]] = userData[k];
					
					//save score for summary file
					userTreeScores[k].push_back(userData[k]);
				}
			}
			
			//get pscores for random trees
			for (int j = 0; j < iters; j++) {
				//create new tree with same num nodes and leaves as users
				randT = new Tree();
				//create random relationships between nodes
				randT->assembleRandomTree();
				//get pscore of random tree
				randomData = pars->getValues(randT);
				
				for(int r = 0; r < numComp; r++) {
					//add trees pscore to map of scores
					it2 = rscoreFreq[r].find(randomData[r]);
					if (it2 != rscoreFreq[r].end()) {//already have that score
						rscoreFreq[r][randomData[r]]++;
					}else{//first time we have seen this score
						rscoreFreq[r][randomData[r]] = 1;
					}
			
					//add randoms score to validscores
					validScores[r][randomData[r]] = randomData[r];
				}
				
				delete randT;
			}
		}else {
			//get pscores for random trees
			for (int j = 0; j < iters; j++) {
				//create new tree with same num nodes and leaves as users
				randT = new Tree();
				//create random relationships between nodes
				randT->assembleRandomTree();
				//get pscore of random tree
				randomData = pars->getValues(randT);
				
				for(int r = 0; r < numComp; r++) {
					//add trees pscore to map of scores
					it2 = rscoreFreq[r].find(randomData[r]);
					if (it2 != rscoreFreq[r].end()) {//already have that score
						rscoreFreq[r][randomData[r]]++;
					}else{//first time we have seen this score
						rscoreFreq[r][randomData[r]] = 1;
					}
			
					//add randoms score to validscores
					validScores[r][randomData[r]] = randomData[r];
				}
				
				delete randT;
			}
		}
		
		float rcumul = 0.0000;
		float ucumul = 0.0000;
		
		for(int a = 0; a < numComp; a++) {
		//this loop fills the cumulative maps and put 0.0000 in the score freq map to make it easier to print.
			for (it = validScores[a].begin(); it != validScores[a].end(); it++) { 
				if (randomtree == "") {
					it2 = uscoreFreq[a].find(it->first);
					//user data has that score 
					if (it2 != uscoreFreq[a].end()) { uscoreFreq[a][it->first] /= T.size(); ucumul+= it2->second;  }
					else { uscoreFreq[a][it->first] = 0.0000; } //no user trees with that score
					//make uCumul map
					uCumul[a][it->first] = ucumul-a;
				}
			
				//make rscoreFreq map and rCumul
				it2 = rscoreFreq[a].find(it->first);
				//get percentage of random trees with that info
				if (it2 != rscoreFreq[a].end()) {  rscoreFreq[a][it->first] /= iters; rcumul+= it2->second;  }
				else { rscoreFreq[a][it->first] = 0.0000; } //no random trees with that score
				rCumul[a][it->first] = rcumul-a;
			}
			
			//find the signifigance of each user trees score when compared to the random trees and save for printing the summary file
			for (int h = 0; h < userTreeScores[a].size(); h++) {
				UScoreSig[a].push_back(rCumul[a][userTreeScores[a][h]]);
			}
		}
		
		printParsimonyFile();
		printUSummaryFile();
		
		//reset globaldata's treemap if you just did random distrib
		if (randomtree != "") { globaldata->gTreemap = savetmap; }
		
		//reset randomTree parameter to ""
		globaldata->setRandomTree("");
		//reset groups parameter
		globaldata->Groups.clear();  globaldata->setGroups("");
		
		return 0;
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ParsimonyCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ParsimonyCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************/
void ParsimonyCommand::printParsimonyFile() {
	try {
		//column headers
		if (randomtree == "") {
			out << "Comb" << '\t' << "Score" << '\t' << "UserFreq" << '\t' << "UserCumul" << '\t' << "RandFreq" << '\t' << "RandCumul" << endl;
		}else {
			out << "Comb" << '\t' << "Score" << '\t' << "RandFreq" << '\t' << "RandCumul" << endl;
		}
		
		//format output
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		for(int a = 0; a < numComp; a++) {
			//print each line
			for (it = validScores[a].begin(); it != validScores[a].end(); it++) { 
				if (randomtree == "") {
					out << setprecision(6)  << groupComb[a] << '\t' << it->first << '\t' << '\t'<< uscoreFreq[a][it->first] << '\t' << uCumul[a][it->first] << '\t' << rscoreFreq[a][it->first] << '\t' << rCumul[a][it->first] << endl; 
				}else{
					out << setprecision(6) << groupComb[a] << '\t' << it->first << '\t' << '\t' << rscoreFreq[a][it->first] << '\t' << rCumul[a][it->first] << endl; 	
				}
			} 
		}
		out.close();
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ParsimonyCommand class Function printParsimonyFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ParsimonyCommand class function printParsimonyFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/***********************************************************/
void ParsimonyCommand::printUSummaryFile() {
	try {
		//column headers
		outSum << "Tree#" << '\t' << "Comb" << '\t'  <<  "ParsScore" << '\t' << '\t' << "ParsSig" <<  endl;
		
		//format output
		outSum.setf(ios::fixed, ios::floatfield); outSum.setf(ios::showpoint);
		
		
		//print each line
		for (int i = 0; i< T.size(); i++) {
			for(int a = 0; a < numComp; a++) {
				outSum << setprecision(6) << i+1 << '\t' << groupComb[a] << '\t' << '\t' << userTreeScores[a][i] << '\t' << UScoreSig[a][i] << endl;
			}
		}
		
		outSum.close();
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ParsimonyCommand class Function printUSummaryFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ParsimonyCommand class function printUSummaryFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************/
void ParsimonyCommand::getUserInput() {
	try {
	
		//create treemap
		tmap = new TreeMap();

		cout << "Please enter the number of groups you would like to analyze: ";
		cin >> numGroups;
			
		int num, count;
		count = 1;
		numEachGroup.resize(numGroups, 0);  
		
		for (int i = 1; i <= numGroups; i++) {
			cout << "Please enter the number of sequences in group " << i <<  ": ";
			cin >> num;
				
			//set tmaps seqsPerGroup
			tmap->seqsPerGroup[toString(i)] = num;
			
			//set tmaps namesOfSeqs
			for (int j = 0; j < num; j++) {
				tmap->namesOfSeqs.push_back(toString(count));
				tmap->treemap[toString(count)].groupname = toString(i);
				count++;
			}
		}
		
		//clears buffer so next command doesn't have error
		string s;	
		getline(cin, s);
		
		//save tmap for later
		globaldata->gTreemap = tmap;
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ParsimonyCommand class Function getUserInput. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ParsimonyCommand class function getUserInput. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
/***********************************************************/

void ParsimonyCommand::setGroups() {
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
						allGroups += tmap->namesOfGroups[i];
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
		groupComb.push_back(allGroups);
		numComp++;
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ParsimonyCommand class Function setGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ParsimonyCommand class function setGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		

}
/*****************************************************************/


