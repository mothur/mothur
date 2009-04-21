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

		util = new SharedUtil();
		util->setGroups(globaldata->Groups, tmap->namesOfGroups, allGroups, numGroups, "unweighted");	//sets the groups the user wants to analyze
		util->getCombos(groupComb, globaldata->Groups, numComp);
		globaldata->setGroups("");
		
		if (numGroups == 1) { numComp++; groupComb.push_back(allGroups); }
				
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
		
		outSum << "Tree#" << '\t' << "Groups" << '\t'  <<  "UWScore" <<'\t' << "UWSig" <<  endl;
		cout << "Tree#" << '\t' << "Groups" << '\t'  <<  "UWScore" << '\t' << "UWSig" <<  endl;
		
		//get pscores for users trees
		for (int i = 0; i < T.size(); i++) {
			counter = 0;
			
			output = new ColumnFile(globaldata->getTreeFile()  + toString(i+1) + ".unweighted");
			
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
			printUWSummaryFile(i);
			
			delete output;
			rscoreFreq.clear(); 
			rCumul.clear();  
			validScores.clear(); 
			utreeScores.clear();  
			UWScoreSig.clear(); 
		}
		
		//reset groups parameter
		globaldata->Groups.clear(); 
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
		vector<string> tags;
		tags.push_back("Score"); tags.push_back("RandFreq"); tags.push_back("RandCumul");
		
		for(int a = 0; a < numComp; a++) {
			output->initFile(groupComb[a], tags);
			//print each line
			for (it = validScores.begin(); it != validScores.end(); it++) { 
				data.push_back(it->first);  data.push_back(rscoreFreq[a][it->first]); data.push_back(rCumul[a][it->first]); 
				output->output(data);
				data.clear();
			} 
			output->resetFile();
		}
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
void UnifracUnweightedCommand::printUWSummaryFile(int i) {
	try {
				
		//format output
		outSum.setf(ios::fixed, ios::floatfield); outSum.setf(ios::showpoint);
			
		//print each line

		for(int a = 0; a < numComp; a++) {
			outSum << i+1 << '\t';
			cout << i+1 << '\t';
			
			if (UWScoreSig[a][0] > (1/(float)iters)) {
				outSum << setprecision(6) << groupComb[a]  << '\t' << utreeScores[a][0] << '\t' << setprecision(globaldata->getIters().length()) << UWScoreSig[a][0] << endl;
				cout << setprecision(6)  << groupComb[a]  << '\t' << utreeScores[a][0] << '\t' << setprecision(globaldata->getIters().length()) << UWScoreSig[a][0] << endl; 
			}else {
				outSum << setprecision(6) << groupComb[a]  << '\t' << utreeScores[a][0] << '\t' << setprecision(globaldata->getIters().length()) << "<" << (1/float(iters)) << endl;
				cout << setprecision(6)  << groupComb[a]  << '\t' << utreeScores[a][0] << '\t' << setprecision(globaldata->getIters().length()) << "<" << (1/float(iters)) << endl; 
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


