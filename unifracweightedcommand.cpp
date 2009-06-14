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
UnifracWeightedCommand::UnifracWeightedCommand(string option) {
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		Groups.clear();
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"groups","iters"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters=parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			if (globaldata->gTree.size() == 0) {//no trees were read
				cout << "You must execute the read.tree command, before you may execute the unifrac.weighted command." << endl; abort = true;  }
										
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				splitAtDash(groups, Groups);
				globaldata->Groups = Groups;
			}
				
			itersString = validParameter.validFile(parameters, "iters", false);			if (itersString == "not found") { itersString = "1000"; }
			convert(itersString, iters); 
		
			
			if (abort == false) {
				T = globaldata->gTree;
				tmap = globaldata->gTreemap;
				sumFile = globaldata->getTreeFile() + ".wsummary";
				openOutputFile(sumFile, outSum);
				
				util = new SharedUtil();
				string s; //to make work with setgroups
				util->setGroups(globaldata->Groups, tmap->namesOfGroups, s, numGroups, "weighted");	//sets the groups the user wants to analyze
				util->getCombos(groupComb, globaldata->Groups, numComp);
				
				weighted = new Weighted(tmap);
				
			}
		}
		
		
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
//**********************************************************************************************************************

void UnifracWeightedCommand::help(){
	try {
		cout << "The unifrac.weighted command can only be executed after a successful read.tree command." << "\n";
		cout << "The unifrac.weighted command parameters are groups and iters.  No parameters are required." << "\n";
		cout << "The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 2 valid groups." << "\n";
		cout << "The group names are separated by dashes.  The iters parameter allows you to specify how many random trees you would like compared to your tree." << "\n";
		cout << "The unifrac.weighted command should be in the following format: unifrac.weighted(groups=yourGroups, iters=yourIters)." << "\n";
		cout << "Example unifrac.weighted(groups=A-B-C, iters=500)." << "\n";
		cout << "The default value for groups is all the groups in your groupfile, and iters is 1000." << "\n";
		cout << "The unifrac.weighted command output two files: .weighted and .wsummary their descriptions are in the manual." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups)." << "\n" << "\n";
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the UnifracWeightedCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the UnifracWeightedCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

/***********************************************************/
int UnifracWeightedCommand::execute() {
	try {
	
		if (abort == true) { return 0; }
		
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
			
			output = new ColumnFile(globaldata->getTreeFile()  + toString(i+1) + ".weighted", itersString);

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
			for (map<float,float>::iterator it = validScores.begin(); it != validScores.end(); it++) { 
				data.push_back(it->first);  data.push_back(rScoreFreq[a][it->first]); data.push_back(rCumul[a][it->first]); 
				output->output(data);
				data.clear();
			} 
			output->resetFile();
		}
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
					outSum << setprecision(6) << i+1 << '\t' << groupComb[j] << '\t' << utreeScores[count] << '\t' << setprecision(itersString.length()) << WScoreSig[count] << endl; 
					cout << setprecision(6) << i+1 << '\t' << groupComb[j] << '\t' << utreeScores[count] << '\t' << setprecision(itersString.length()) << WScoreSig[count] << endl; 
				}else{
					outSum << setprecision(6) << i+1 << '\t' << groupComb[j] << '\t' << utreeScores[count] << '\t' << setprecision(itersString.length()) << "<" << (1/float(iters)) << endl; 
					cout << setprecision(6) << i+1 << '\t' << groupComb[j] << '\t' << utreeScores[count] << '\t' << setprecision(itersString.length()) << "<" << (1/float(iters)) << endl; 
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
				map<float,float>::iterator it = rScoreFreq[f].find(rScores[f][i]);
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
			for (map<float,float>::iterator it = validScores.begin(); it != validScores.end(); it++) {
				//make rscoreFreq map and rCumul
				map<float,float>::iterator it2 = rScoreFreq[a].find(it->first);
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

/***********************************************************/





