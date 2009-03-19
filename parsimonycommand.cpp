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
			parsFileout = globaldata->getTreeFile() + "temp" + ".parsimony";
			sumFile = globaldata->getTreeFile() + ".psummary";
			openOutputFile(sumFile, outSum);
		}else { //user wants random distribution
			savetmap = globaldata->gTreemap;
			getUserInput();
			parsFile = randomtree;
			parsFileout = globaldata->getTreeFile() + "temp";
		}
		
		//set users groups to analyze
		setGroups();
		convert(globaldata->getIters(), iters);  //how many random trees to generate
		pars = new Parsimony(tmap);
		counter = 0;

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
		userData.resize(numComp,0);  //data = AB, AC, BC, ABC.
		randomData.resize(numComp,0);  //data = AB, AC, BC, ABC.
		rscoreFreq.resize(numComp);  
		uscoreFreq.resize(numComp);  
		rCumul.resize(numComp);  
		uCumul.resize(numComp);  
		userTreeScores.resize(numComp);  
		UScoreSig.resize(numComp); 
				
		if (randomtree == "") {
			//get pscores for users trees
			for (int i = 0; i < T.size(); i++) {
				userData = pars->getValues(T[i]);  //data = AB, AC, BC, ABC.
				
				//output scores for each combination
				for(int k = 0; k < numComp; k++) {
					//update uscoreFreq
					it = uscoreFreq[k].find(userData[k]);
					if (it == uscoreFreq[k].end()) {//new score
						uscoreFreq[k][userData[k]] = 1;
					}else{ uscoreFreq[k][userData[k]]++; }
					
					//add users score to valid scores
					validScores[userData[k]] = userData[k];
					
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
					validScores[randomData[r]] = randomData[r];
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
					validScores[randomData[r]] = randomData[r];
				}
				
				delete randT;
			}
		}
		
		for(int a = 0; a < numComp; a++) {
			float rcumul = 0.0000;
			float ucumul = 0.0000;
			//this loop fills the cumulative maps and put 0.0000 in the score freq map to make it easier to print.
			for (it = validScores.begin(); it != validScores.end(); it++) { 
				if (randomtree == "") {
					it2 = uscoreFreq[a].find(it->first);
					//user data has that score 
					if (it2 != uscoreFreq[a].end()) { uscoreFreq[a][it->first] /= T.size(); ucumul+= it2->second;  }
					else { uscoreFreq[a][it->first] = 0.0000; } //no user trees with that score
					//make uCumul map
					uCumul[a][it->first] = ucumul;
				}
			
				//make rscoreFreq map and rCumul
				it2 = rscoreFreq[a].find(it->first);
				//get percentage of random trees with that info
				if (it2 != rscoreFreq[a].end()) {  rscoreFreq[a][it->first] /= iters; rcumul+= it2->second;  }
				else { rscoreFreq[a][it->first] = 0.0000; } //no random trees with that score
				rCumul[a][it->first] = rcumul;
			}
			
			//find the signifigance of each user trees score when compared to the random trees and save for printing the summary file
			for (int h = 0; h < userTreeScores[a].size(); h++) {
				UScoreSig[a].push_back(rCumul[a][userTreeScores[a][h]]);
			}
		}
		
		printParsimonyFile();
		if (randomtree == "") { printUSummaryFile(); }
		
		//reset globaldata's treemap if you just did random distrib
		if (randomtree != "") {
			//memory leak prevention
			if (globaldata->gTreemap != NULL) { delete globaldata->gTreemap;  }
			globaldata->gTreemap = savetmap;
		}
		
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
		vector<double> data;
		
		//format output
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);

		for(int a = 0; a < numComp; a++) {
			initFile(groupComb[a]);
			//print each line
			for (it = validScores.begin(); it != validScores.end(); it++) { 
				if (randomtree == "") {
					data.push_back(it->first);  data.push_back(uscoreFreq[a][it->first]); data.push_back(uCumul[a][it->first]); data.push_back(rscoreFreq[a][it->first]); data.push_back(rCumul[a][it->first]); 
				}else{
					data.push_back(it->first);  data.push_back(rscoreFreq[a][it->first]); data.push_back(rCumul[a][it->first]); 
				}
				output(data);
				data.clear();
			} 
			resetFile();
		}
		
		out.close();
		inFile.close();
		remove(parsFileout.c_str());
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
		outSum << "Tree#" << '\t' << "Groups" << '\t'  <<  "ParsScore" << '\t' << "ParsSig" <<  endl;
		cout << "Tree#" << '\t' << "Groups" << '\t'  <<  "ParsScore" << '\t' << "ParsSig" <<  endl;
		
		//format output
		outSum.setf(ios::fixed, ios::floatfield); outSum.setf(ios::showpoint);
		
		
		//print each line
		for (int i = 0; i< T.size(); i++) {
			for(int a = 0; a < numComp; a++) {
				if (UScoreSig[a][i] > (1/(float)iters)) {
					outSum << setprecision(6) << i+1 << '\t' << groupComb[a]  << '\t' << userTreeScores[a][i] << setprecision(globaldata->getIters().length()) << '\t' << UScoreSig[a][i] << endl;
					cout << setprecision(6) << i+1 << '\t' << groupComb[a]  << '\t' << userTreeScores[a][i] << setprecision(globaldata->getIters().length()) << '\t' << UScoreSig[a][i] << endl;
				}else {
					outSum << setprecision(6) << i+1 << '\t' << groupComb[a] << '\t' << userTreeScores[a][i] << setprecision(globaldata->getIters().length())  << '\t' << "<" << (1/float(iters)) << endl;
					cout << setprecision(6) << i+1 << '\t' << groupComb[a] << '\t' << userTreeScores[a][i] << setprecision(globaldata->getIters().length()) << '\t' << "<" << (1/float(iters)) << endl;
				}
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
			tmap->namesOfGroups.push_back(toString(i));
			
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
		//memory leak prevention
		if (globaldata->gTreemap != NULL) { delete globaldata->gTreemap;  }
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
				groupComb.push_back(globaldata->Groups[r]+ "-" +globaldata->Groups[l]);
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
		cout << "Standard Error: " << e.what() << " has occurred in the ParsimonyCommand class Function setGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ParsimonyCommand class function setGroups. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}		

}
/*****************************************************************/

void ParsimonyCommand::initFile(string label){
	try {
		if(counter != 0){
			openOutputFile(parsFileout, out);
			openInputFile(parsFile, inFile);

			string inputBuffer;
			getline(inFile, inputBuffer);
			
			if (randomtree == "") {
				out <<  inputBuffer << '\t' << label + "Score" << '\t' << label + "UserFreq" << '\t' << label + "UserCumul" << '\t' << label + "RandFreq" << '\t' << label + "RandCumul" << endl;
			}else {
				out <<  inputBuffer << '\t' << "Score" << '\t' << "RandFreq" << '\t' << "RandCumul" << endl;
			}
		}else{
			openOutputFile(parsFileout, out);
			//column headers
			if (randomtree == "") {
				out << label + "Score" << '\t' << label + "UserFreq" << '\t' << label + "UserCumul" << '\t' << label + "RandFreq" << '\t' << label + "RandCumul" << endl;
			}else {
				out << "Score" << '\t' << "RandFreq" << '\t' << "RandCumul" << endl;
			}
		}

		out.setf(ios::fixed, ios::floatfield);
		out.setf(ios::showpoint);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ParsimonyCommand class Function initFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ParsimonyCommand class function initFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

void ParsimonyCommand::output(vector<double> data){
	try {
		if(counter != 0){		
			string inputBuffer;
			getline(inFile, inputBuffer);
		
			if (randomtree == "") {
				out << inputBuffer << '\t' << setprecision(6) << data[0] << setprecision(globaldata->getIters().length())  << '\t' << data[1] << '\t' << data[2] << '\t' << data[3] << '\t' << data[4] << endl;
			}else{
				out << inputBuffer << '\t' << setprecision(6) << data[0] << setprecision(globaldata->getIters().length())  << '\t' << data[1] << '\t' << data[2] << endl;
			}
		}
		else{
			if (randomtree == "") {
				out << setprecision(6) << data[0] << setprecision(globaldata->getIters().length())  << '\t' << data[1] << '\t' << data[2] << '\t' << data[3] << '\t' << data[4] << endl;
			}else{
				out << setprecision(6) << data[0] << setprecision(globaldata->getIters().length())  << '\t' << data[1] << '\t' << data[2] << endl;
			}
		}

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ParsimonyCommand class Function output. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ParsimonyCommand class function output. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************/

void ParsimonyCommand::resetFile(){
	try {
		if(counter != 0){
			out.close();
			inFile.close();
		}
		else{
			out.close();
		}
		counter = 1;
		
		remove(parsFile.c_str());
		rename(parsFileout.c_str(), parsFile.c_str());
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ParsimonyCommand class Function resetFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ParsimonyCommand class function resetFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}


