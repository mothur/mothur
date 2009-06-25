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
UnifracUnweightedCommand::UnifracUnweightedCommand(string option) {
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
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			if (globaldata->gTree.size() == 0) {//no trees were read
				mothurOut("You must execute the read.tree command, before you may execute the unifrac.unweighted command."); mothurOutEndLine(); abort = true;  }
										
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
				sumFile = globaldata->getTreeFile() + ".uwsummary";
				openOutputFile(sumFile, outSum);
				
				util = new SharedUtil();
				util->setGroups(globaldata->Groups, tmap->namesOfGroups, allGroups, numGroups, "unweighted");	//sets the groups the user wants to analyze
				util->getCombos(groupComb, globaldata->Groups, numComp);
				
				if (numGroups == 1) { numComp++; groupComb.push_back(allGroups); }
				
				unweighted = new Unweighted(tmap);
				
			}
			
		}
		
	}
	catch(exception& e) {
		errorOut(e, "UnifracUnweightedCommand", "UnifracUnweightedCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void UnifracUnweightedCommand::help(){
	try {
		mothurOut("The unifrac.unweighted command can only be executed after a successful read.tree command.\n");
		mothurOut("The unifrac.unweighted command parameters are groups and iters.  No parameters are required.\n");
		mothurOut("The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 1 valid group.\n");
		mothurOut("The group names are separated by dashes.  The iters parameter allows you to specify how many random trees you would like compared to your tree.\n");
		mothurOut("The unifrac.unweighted command should be in the following format: unifrac.unweighted(groups=yourGroups, iters=yourIters).\n");
		mothurOut("Example unifrac.unweighted(groups=A-B-C, iters=500).\n");
		mothurOut("The default value for groups is all the groups in your groupfile, and iters is 1000.\n");
		mothurOut("The unifrac.unweighted command output two files: .unweighted and .uwsummary their descriptions are in the manual.\n");
		mothurOut("Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n\n");
	}
	catch(exception& e) {
		errorOut(e, "UnifracUnweightedCommand", "help");
		exit(1);
	}
}


/***********************************************************/
int UnifracUnweightedCommand::execute() {
	try {
		
		if (abort == true) { return 0; }
		
		userData.resize(numComp,0);  //data[0] = unweightedscore 
		randomData.resize(numComp,0); //data[0] = unweightedscore
		//create new tree with same num nodes and leaves as users
		
		outSum << "Tree#" << '\t' << "Groups" << '\t'  <<  "UWScore" <<'\t' << "UWSig" <<  endl;
		mothurOut("Tree#\tGroups\tUWScore\tUWSig"); mothurOutEndLine();
		
		//get pscores for users trees
		for (int i = 0; i < T.size(); i++) {
			counter = 0;
			
			output = new ColumnFile(globaldata->getTreeFile()  + toString(i+1) + ".unweighted", itersString);
			
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
					map<float,float>::iterator it = rscoreFreq[k].find(randomData[k]);
					if (it != rscoreFreq[k].end()) {//already have that score
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
				for (map<float,float>::iterator it = validScores.begin(); it != validScores.end(); it++) { 
					//make rscoreFreq map and rCumul
					map<float,float>::iterator it2 = rscoreFreq[a].find(it->first);
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
		errorOut(e, "UnifracUnweightedCommand", "execute");
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
			for (map<float,float>::iterator it = validScores.begin(); it != validScores.end(); it++) { 
				data.push_back(it->first);  data.push_back(rscoreFreq[a][it->first]); data.push_back(rCumul[a][it->first]); 
				output->output(data);
				data.clear();
			} 
			output->resetFile();
		}
	}
	catch(exception& e) {
		errorOut(e, "UnifracUnweightedCommand", "printUnweightedFile");
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
			mothurOut(toString(i+1) + "\t");
			
			if (UWScoreSig[a][0] > (1/(float)iters)) {
				outSum << setprecision(6) << groupComb[a]  << '\t' << utreeScores[a][0] << '\t' << setprecision(itersString.length()) << UWScoreSig[a][0] << endl;
				cout << setprecision(6)  << groupComb[a]  << '\t' << utreeScores[a][0] << '\t' << setprecision(itersString.length()) << UWScoreSig[a][0] << endl; 
				mothurOutJustToLog(groupComb[a]  + "\t" + toString(utreeScores[a][0])  + "\t" + toString(UWScoreSig[a][0])); mothurOutEndLine(); 
			}else {
				outSum << setprecision(6) << groupComb[a]  << '\t' << utreeScores[a][0] << '\t' << setprecision(itersString.length()) << "<" << (1/float(iters)) << endl;
				cout << setprecision(6)  << groupComb[a]  << '\t' << utreeScores[a][0] << '\t' << setprecision(itersString.length()) << "<" << (1/float(iters)) << endl; 
				mothurOutJustToLog(groupComb[a]  + "\t" + toString(utreeScores[a][0])  + "\t<" + toString((1/float(iters)))); mothurOutEndLine();
			}
		}
		
	}
	catch(exception& e) {
		errorOut(e, "UnifracUnweightedCommand", "printUWSummaryFile");
		exit(1);
	}
}

/***********************************************************/


