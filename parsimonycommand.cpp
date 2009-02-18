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
			distFile = globaldata->getTreeFile() + ".pdistrib";
			openOutputFile(distFile, outDist);
			
			//if the user has not entered specific groups to analyze then do them all
			if (globaldata->Groups.size() != 0) {
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
					cout << "When using the groups parameter you must have at least 1 valid group. I will run the command using all the groups in your groupfile." << endl; 
					for (int i = 0; i < tmap->namesOfGroups.size(); i++) {
						globaldata->Groups.push_back(tmap->namesOfGroups[i]);
					}
				}		
			}else {
				for (int i = 0; i < tmap->namesOfGroups.size(); i++) {
					globaldata->Groups.push_back(tmap->namesOfGroups[i]);
				}
			}

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
		
		//get pscore for users tree
		userData.resize(1,0);  //data[0] = pscore.
		randomData.resize(1,0);  //data[0] = pscore.
		
		//format output
		outDist.setf(ios::fixed, ios::floatfield); outDist.setf(ios::showpoint);
		outDist << "RandomTree#" << '\t' << "ParsScore" << endl;
		
		if (randomtree == "") {
			//get pscores for users trees
			for (int i = 0; i < T.size(); i++) {
				cout << "Processing tree " << i+1 << endl;
				userData = pars->getValues(T[i]);  //userData[0] = pscore
				cout << "Tree " << i+1 << " parsimony score = " << userData[0] << endl;
				//update uscoreFreq
				it = uscoreFreq.find(userData[0]);
				if (it == uscoreFreq.end()) {//new score
					uscoreFreq[userData[0]] = 1;
				}else{ uscoreFreq[userData[0]]++; }
			
				//add users score to valid scores
				validScores[userData[0]] = userData[0];
				
				//save score for summary file
				userTreeScores.push_back(userData[0]);
			
			}
			
			//get pscores for random trees
			for (int j = 0; j < iters; j++) {
				//create new tree with same num nodes and leaves as users
				randT = new Tree();
				//create random relationships between nodes
				randT->assembleRandomTree();
				//get pscore of random tree
				randomData = pars->getValues(randT);
			
				//add trees pscore to map of scores
				it2 = rscoreFreq.find(randomData[0]);
				if (it2 != rscoreFreq.end()) {//already have that score
					rscoreFreq[randomData[0]]++;
				}else{//first time we have seen this score
					rscoreFreq[randomData[0]] = 1;
				}
			
				//add randoms score to validscores
				validScores[randomData[0]] = randomData[0];
				
				//output info to pdistrib file
				outDist << j+1 << '\t'<< '\t' << randomData[0] << endl;
					
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
			
				//add trees pscore to map of scores
				it2 = rscoreFreq.find(randomData[0]);
				if (it2 != rscoreFreq.end()) {//already have that score
					rscoreFreq[randomData[0]]++;
				}else{//first time we have seen this score
					rscoreFreq[randomData[0]] = 1;
				}
			
				//add randoms score to validscores
				validScores[randomData[0]] = randomData[0];
					
				delete randT;
			}
		}
		
		float rcumul = 0.0000;
		float ucumul = 0.0000;
		
		//this loop fills the cumulative maps and put 0.0000 in the score freq map to make it easier to print.
		for (it = validScores.begin(); it != validScores.end(); it++) { 
			if (randomtree == "") {
				it2 = uscoreFreq.find(it->first);
				//user data has that score 
				if (it2 != uscoreFreq.end()) { uscoreFreq[it->first] /= T.size(); ucumul+= it2->second;  }
				else { uscoreFreq[it->first] = 0.0000; } //no user trees with that score
				//make uCumul map
				uCumul[it->first] = ucumul;
			}
			
			//make rscoreFreq map and rCumul
			it2 = rscoreFreq.find(it->first);
			//get percentage of random trees with that info
			if (it2 != rscoreFreq.end()) {  rscoreFreq[it->first] /= iters; rcumul+= it2->second;  }
			else { rscoreFreq[it->first] = 0.0000; } //no random trees with that score
			rCumul[it->first] = rcumul;
		}
		
		//find the signifigance of each user trees score when compared to the random trees and save for printing the summary file
		for (int h = 0; h < userTreeScores.size(); h++) {
			UScoreSig.push_back(rCumul[userTreeScores[h]]);
		}

		printParsimonyFile();
		printUSummaryFile();
		
		//reset globaldata's treemap if you just did random distrib
		if (randomtree != "") { globaldata->gTreemap = savetmap; }
		
		//reset randomTree parameter to ""
		globaldata->setRandomTree("");
		//reset groups parameter
		globaldata->Groups.clear();
		
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
			out << "Score" << '\t' << "UserFreq" << '\t' << "UserCumul" << '\t' << "RandFreq" << '\t' << "RandCumul" << endl;
		}else {
			out << "Score" << '\t' << "RandFreq" << '\t' << "RandCumul" << endl;
		}
		
		//format output
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		//print each line
		for (it = validScores.begin(); it != validScores.end(); it++) { 
			if (randomtree == "") {
				out << setprecision(6) << it->first << '\t' << '\t' << uscoreFreq[it->first] << '\t' << uCumul[it->first] << '\t' << rscoreFreq[it->first] << '\t' << rCumul[it->first] << endl; 
			}else{
				out << setprecision(6) << it->first << '\t' << '\t' << rscoreFreq[it->first] << '\t' << rCumul[it->first] << endl; 	
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
		outSum << "Tree#" << '\t'  <<  "ParsScore" << '\t' << '\t' << "ParsSig" <<  endl;
		
		//format output
		outSum.setf(ios::fixed, ios::floatfield); outSum.setf(ios::showpoint);
		
		//print each line
		for (int i = 0; i< T.size(); i++) {
			outSum << setprecision(6) << i+1 << '\t' << '\t' << userTreeScores[i] << '\t' << UScoreSig[i] << endl; 
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

