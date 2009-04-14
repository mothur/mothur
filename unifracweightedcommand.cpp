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
		sumFile = globaldata->getTreeFile() + ".wsummary";
		openOutputFile(sumFile, outSum);
				
		util = new SharedUtil();
		string s; //to make work with setgroups
		util->setGroups(globaldata->Groups, tmap->namesOfGroups, s, numGroups, "weighted");	//sets the groups the user wants to analyze
		util->getCombos(groupComb, globaldata->Groups, numComp);
		globaldata->setGroups("");
				
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
		Progress* reading;
		reading = new Progress("Comparing to random:", iters);
		
		//get weighted for users tree
		userData.resize(numComp,0);  //data[0] = weightedscore AB, data[1] = weightedscore AC...
		randomData.resize(numComp,0); //data[0] = weightedscore AB, data[1] = weightedscore AC...
				
		//create new tree with same num nodes and leaves as users
		randT = new Tree();
		
		//get weighted scores for users trees
		for (int i = 0; i < T.size(); i++) {
			counter = 0;
			rScores.resize(numComp);  //data[0] = weightedscore AB, data[1] = weightedscore AC...
			uScores.resize(numComp);  //data[0] = weightedscore AB, data[1] = weightedscore AC...
			//weightedFile = globaldata->getTreeFile()  + toString(i+1) + ".weighted";
			//weightedFileout = globaldata->getTreeFile() + "temp." + toString(i+1) + ".weighted";
			output = new ThreeColumnFile2(globaldata->getTreeFile()  + toString(i+1) + ".weighted");

			userData = weighted->getValues(T[i]);  //userData[0] = weightedscore
			
			//save users score
			for (int s=0; s<numComp; s++) {
				//add users score to vector of user scores
				uScores[s].push_back(userData[s]);
				
				//save users tree score for summary file
				utreeScores.push_back(userData[s]);
			}
			
			//get scores for random trees
			for (int j = 0; j < iters; j++) {
				int count = 0;
				for (int r=0; r<numGroups; r++) { 
					for (int l = r+1; l < numGroups; l++) {
						//copy T[i]'s info.
						randT->getCopy(T[i]);
						 
						//create a random tree with same topology as T[i], but different labels
						randT->assembleRandomUnifracTree(globaldata->Groups[r], globaldata->Groups[l]);
						//get wscore of random tree
						randomData = weighted->getValues(randT, globaldata->Groups[r], globaldata->Groups[l]);
						
						//save scores
						rScores[count].push_back(randomData[0]);
						count++;
					}
				}
				
				//update progress bar
				reading->update(j);

			}

			//removeValidScoresDuplicates(); 
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
			
			//out << "Tree# " << i << endl;
			calculateFreqsCumuls();
			printWeightedFile();
			
			delete output;
			
			//clear data
			rScores.clear();
			uScores.clear();
			validScores.clear();
		}
		
		//finish progress bar
		reading->finish();
		delete reading;
		
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
/***********************************************************/
void UnifracWeightedCommand::printWeightedFile() {
	try {
		vector<double> data;
		vector<string> tags;
		tags.push_back("Score"); tags.push_back("RandFreq"); tags.push_back("RandCumul");
		
		for(int a = 0; a < numComp; a++) {
			output->initFile(groupComb[a], tags);
			//print each line
			for (it = validScores.begin(); it != validScores.end(); it++) { 
				data.push_back(it->first);  data.push_back(rScoreFreq[a][it->first]); data.push_back(rCumul[a][it->first]); 
				output->output(data);
				data.clear();
			} 
			output->resetFile();
		}
		
		//out.close();
		//inFile.close();
		//remove(weightedFileout.c_str());
		
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
		outSum << "Tree#" << '\t' << "Groups" << '\t' << "WScore" << '\t' << "WSig" <<  endl;
		cout << "Tree#" << '\t' << "Groups" << '\t' << "WScore" << '\t' << "WSig" <<  endl;
		
		//format output
		outSum.setf(ios::fixed, ios::floatfield); outSum.setf(ios::showpoint);
		
		//print each line
		int count = 0;
		for (int i = 0; i < T.size(); i++) { 
			for (int j = 0; j < numComp; j++) {
				if (WScoreSig[count] > (1/(float)iters)) {
					outSum << setprecision(6) << i+1 << '\t' << groupComb[j] << '\t' << utreeScores[count] << '\t' << setprecision(globaldata->getIters().length()) << WScoreSig[count] << endl; 
					cout << setprecision(6) << i+1 << '\t' << groupComb[j] << '\t' << utreeScores[count] << '\t' << setprecision(globaldata->getIters().length()) << WScoreSig[count] << endl; 
				}else{
					outSum << setprecision(6) << i+1 << '\t' << groupComb[j] << '\t' << utreeScores[count] << '\t' << setprecision(globaldata->getIters().length()) << "<" << (1/float(iters)) << endl; 
					cout << setprecision(6) << i+1 << '\t' << groupComb[j] << '\t' << utreeScores[count] << '\t' << setprecision(globaldata->getIters().length()) << "<" << (1/float(iters)) << endl; 
				}
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

void UnifracWeightedCommand::calculateFreqsCumuls() {
	try {
		//clear out old tree values
		rScoreFreq.clear();
		rScoreFreq.resize(numComp);
		rCumul.clear();
		rCumul.resize(numComp);
		validScores.clear();
	
		//calculate frequency
		for (int f = 0; f < numComp; f++) {
			for (int i = 0; i < rScores[f].size(); i++) { //looks like 0,0,1,1,1,2,4,7...  you want to make a map that say rScoreFreq[0] = 2, rScoreFreq[1] = 3...
				validScores[rScores[f][i]] = rScores[f][i];
				it = rScoreFreq[f].find(rScores[f][i]);
				if (it != rScoreFreq[f].end()) {
					rScoreFreq[f][rScores[f][i]]++;
				}else{
					rScoreFreq[f][rScores[f][i]] = 1;
				}
			}
		}
		
		//calculate rcumul
		for(int a = 0; a < numComp; a++) {
			float rcumul = 1.0000;
			//this loop fills the cumulative maps and put 0.0000 in the score freq map to make it easier to print.
			for (it = validScores.begin(); it != validScores.end(); it++) {
				//make rscoreFreq map and rCumul
				it2 = rScoreFreq[a].find(it->first);
				rCumul[a][it->first] = rcumul;
				//get percentage of random trees with that info
				if (it2 != rScoreFreq[a].end()) {  rScoreFreq[a][it->first] /= iters; rcumul-= it2->second;  }
				else { rScoreFreq[a][it->first] = 0.0000; } //no random trees with that score
			}
		}

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UnifracWeightedCommand class Function calculateFreqsCums. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UnifracWeightedCommand class function calculateFreqsCums. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}

/*****************************************************************

void UnifracWeightedCommand::initFile(string label){
	try {
		if(counter != 0){
			openOutputFile(weightedFileout, out);
			openInputFile(weightedFile, inFile);

			string inputBuffer;
			getline(inFile, inputBuffer);
		
			out	<<  inputBuffer << '\t' << label + "RandFreq" << '\t' << label + "RandCumul" << endl;		
		}else{
			openOutputFile(weightedFileout, out);
			out	<< label + "Score" << '\t' << label + "RandFreq" << '\t' << label + "RandCumul" << endl;
		}

		out.setf(ios::fixed, ios::floatfield);
		out.setf(ios::showpoint);
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UnifracWeightedCommand class Function initFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UnifracWeightedCommand class function initFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

/***********************************************************************

void UnifracWeightedCommand::output(vector<double> data){
	try {
		if(counter != 0){		
			string inputBuffer;
			getline(inFile, inputBuffer);

			out << inputBuffer << '\t' << setprecision(6) << data[0] << setprecision(globaldata->getIters().length())  << '\t' << data[1] << '\t' << data[2] << endl;
		}
		else{
			out << setprecision(6) << data[0] << setprecision(globaldata->getIters().length())  << '\t' << data[1] << '\t' << data[2] << endl;

		}
		
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UnifracWeightedCommand class Function output. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UnifracWeightedCommand class function output. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
};

/***********************************************************************

void UnifracWeightedCommand::resetFile(){
	try {
		if(counter != 0){
			out.close();
			inFile.close();
		}
		else{
			out.close();
		}
		counter = 1;
		
		remove(weightedFile.c_str());
		rename(weightedFileout.c_str(), weightedFile.c_str());
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UnifracWeightedCommand class Function resetFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UnifracWeightedCommand class function resetFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	
	
			
}

*/




