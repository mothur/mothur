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
		sumFile = globaldata->getTreeFile() + ".uwsummary";
		openOutputFile(sumFile, outSum);
		distFile = globaldata->getTreeFile() + ".uwdistrib";
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
			}		
		}

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
		
		//get unweighted for users tree
		userData.resize(1,0);  //data[0] = unweightedscore 
		randomData.resize(1,0); //data[0] = unweightedscore
		
		//format output
		outDist.setf(ios::fixed, ios::floatfield); outDist.setf(ios::showpoint);
		outDist << "Tree#" << '\t' << "Iter" << '\t' << "UWScore" << endl;
		
		//create new tree with same num nodes and leaves as users
		randT = new Tree();
 			
		//get pscores for users trees
		for (int i = 0; i < T.size(); i++) {
			cout << "Processing tree " << i+1 << endl;
			userData = unweighted->getValues(T[i]);  //userData[0] = unweightedscore
			
			//update uscoreFreq
			it = uscoreFreq.find(userData[0]);
			if (it == uscoreFreq.end()) {//new score
				uscoreFreq[userData[0]] = 1;
			}else{ uscoreFreq[userData[0]]++; }
			
			//add users score to valid scores
			validScores[userData[0]] = userData[0];
			
			//saves users score
			utreeScores.push_back(userData[0]);
			
			//copy T[i]'s info.
			randT->getCopy(T[i]); 
			
			//get unweighted scores for random trees
			for (int j = 0; j < iters; j++) {
				//create a random tree with same topology as T[i], but different labels
				randT->assembleRandomUnifracTree();
				//get pscore of random tree
				randomData = unweighted->getValues(randT);
			
				//add trees unweighted score to map of scores
				it2 = rscoreFreq.find(randomData[0]);
				if (it2 != rscoreFreq.end()) {//already have that score
					rscoreFreq[randomData[0]]++;
				}else{//first time we have seen this score
					rscoreFreq[randomData[0]] = 1;
				}
				
				//add randoms score to validscores
				validScores[randomData[0]] = randomData[0];
				
				//output info to uwdistrib file
				outDist << i+1 << '\t' << '\t'<< j+1 << '\t' << '\t' << randomData[0] << endl;
			}
			
			saveRandomScores(); //save all random scores for unweighted file
			
			//find the signifigance of the score
			float rcumul = 0.0000;
			for (it = rscoreFreq.begin(); it != rscoreFreq.end(); it++) { 
				//get percentage of random trees with that info
				rscoreFreq[it->first] /= iters; 
				rcumul+= it->second;  
				rCumul[it->first] = rcumul;
			}
			
			//save the signifigance of the users score for printing later
			UWScoreSig.push_back(rCumul[userData[0]]);
			
			
			//clear random data
			rscoreFreq.clear();  //you clear this because in the summary file you want the unweighted signifinance to be relative to these 1000 trees.
			rCumul.clear();
		}
		
		float ucumul = 0.0000;
		float rcumul = 0.0000;
		//this loop fills the cumulative maps and put 0.0000 in the score freq map to make it easier to print.
		for (it = validScores.begin(); it != validScores.end(); it++) { 
			it2 = uscoreFreq.find(it->first);
			//user data has that score 
			if (it2 != uscoreFreq.end()) { uscoreFreq[it->first] /= T.size(); ucumul+= it2->second;  }
			else { uscoreFreq[it->first] = 0.0000; } //no user trees with that score
			//make uCumul map
			uCumul[it->first] = ucumul;
			
			//make rscoreFreq map and rCumul
			it2 = totalrscoreFreq.find(it->first);
			//get percentage of random trees with that info
			if (it2 != totalrscoreFreq.end()) {  totalrscoreFreq[it->first] /= (iters*T.size()); rcumul+= it2->second;  }
			else { totalrscoreFreq[it->first] = 0.0000; } //no random trees with that score
			rCumul[it->first] = rcumul;
		}
		
		printUnweightedFile();
		printUWSummaryFile();
		
		//reset randomTree parameter to 0
		globaldata->setRandomTree("0");
		
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
		//column headers
		
		out << "Score" << '\t' << "UserFreq" << '\t' << "UserCumul" << '\t' << "RandFreq" << '\t' << "RandCumul" << endl;
				
		//format output
		out.setf(ios::fixed, ios::floatfield); out.setf(ios::showpoint);
		
		//print each line
		for (it = validScores.begin(); it != validScores.end(); it++) { 
			out << setprecision(6) << it->first << '\t' << '\t' << uscoreFreq[it->first] << '\t' << uCumul[it->first] << '\t' << totalrscoreFreq[it->first] << '\t' << rCumul[it->first] << endl; 
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
		//column headers
		outSum << "Tree#" << '\t'  <<  "UWScore" << '\t' << '\t' << "UWSig" <<  endl;
		
		//format output
		outSum.setf(ios::fixed, ios::floatfield); outSum.setf(ios::showpoint);
		
		//print each line
		for (int i = 0; i< T.size(); i++) {
			outSum << setprecision(6) << i+1 << '\t' << '\t' << utreeScores[i] << '\t' << UWScoreSig[i] << endl; 
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
void UnifracUnweightedCommand::saveRandomScores() {
	try {
		for (it = rscoreFreq.begin(); it != rscoreFreq.end(); it++) { 
			//does this score already exist in the total map
			it2 = totalrscoreFreq.find(it->first);
			//if yes then add them
			if (it2 != totalrscoreFreq.end()) { 
				totalrscoreFreq[it->first] += rscoreFreq[it->first];
			}else{ //its a new score
				totalrscoreFreq[it->first] = rscoreFreq[it->first];
			}
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UnifracUnweightedCommand class Function saveRandomScores. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UnifracUnweightedCommand class function saveRandomScores. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************/