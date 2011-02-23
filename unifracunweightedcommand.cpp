/*
 *  unifracunweightedcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 2/9/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "unifracunweightedcommand.h"

//**********************************************************************************************************************
vector<string> UnifracUnweightedCommand::getValidParameters(){	
	try {
		string Array[] =  {"groups","iters","distance","random", "processors","outputdir","inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
UnifracUnweightedCommand::UnifracUnweightedCommand(){	
	try {
		globaldata = GlobalData::getInstance();
		abort = true; calledHelp = true; 
		vector<string> tempOutNames;
		outputTypes["unweighted"] = tempOutNames;
		outputTypes["uwsummary"] = tempOutNames;
		outputTypes["phylip"] = tempOutNames;
		outputTypes["column"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "UnifracUnweightedCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> UnifracUnweightedCommand::getRequiredParameters(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> UnifracUnweightedCommand::getRequiredFiles(){	
	try {
		string Array[] =  {"tree","group"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));

		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "getRequiredFiles");
		exit(1);
	}
}
/***********************************************************/
UnifracUnweightedCommand::UnifracUnweightedCommand(string option)  {
	try {
		globaldata = GlobalData::getInstance();
		abort = false; calledHelp = false;   
		Groups.clear();
			
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"groups","iters","distance","random", "processors","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["unweighted"] = tempOutNames;
			outputTypes["uwsummary"] = tempOutNames;
			outputTypes["phylip"] = tempOutNames;
			outputTypes["column"] = tempOutNames;
			
			if (globaldata->gTree.size() == 0) {//no trees were read
				m->mothurOut("You must execute the read.tree command, before you may execute the unifrac.unweighted command."); m->mothurOutEndLine(); abort = true;  }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(globaldata->inputFileName); //if user entered a file with a path then preserve it	
			}
							
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				m->splitAtDash(groups, Groups);
				globaldata->Groups = Groups;
			}
				
			itersString = validParameter.validFile(parameters, "iters", false);				if (itersString == "not found") { itersString = "1000"; }
			convert(itersString, iters); 
			
			string temp = validParameter.validFile(parameters, "distance", false);			
			if (temp == "not found") { phylip = false; outputForm = ""; }
			else{
				if ((temp == "lt") || (temp == "column") || (temp == "square")) {  phylip = true;  outputForm = temp; }
				else { m->mothurOut("Options for distance are: lt, square, or column. Using lt."); m->mothurOutEndLine(); phylip = true; outputForm = "lt"; }
			}
			
			temp = validParameter.validFile(parameters, "random", false);					if (temp == "not found") { temp = "f"; }
			random = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = "1";				}
			convert(temp, processors); 
			
			if (!random) {  iters = 0;  } //turn off random calcs
			
			//if user selects distance = true and no groups it won't calc the pairwise
			if ((phylip) && (Groups.size() == 0)) {
				groups = "all";
				m->splitAtDash(groups, Groups);
				globaldata->Groups = Groups;
			}
		
			if (abort == false) {
				T = globaldata->gTree;
				tmap = globaldata->gTreemap;
				sumFile = outputDir + m->getSimpleName(globaldata->getTreeFile()) + ".uwsummary";
				outputNames.push_back(sumFile); outputTypes["uwsummary"].push_back(sumFile);
				m->openOutputFile(sumFile, outSum);
				
				util = new SharedUtil();
				util->setGroups(globaldata->Groups, tmap->namesOfGroups, allGroups, numGroups, "unweighted");	//sets the groups the user wants to analyze
				util->getCombos(groupComb, globaldata->Groups, numComp);
				
				if (numGroups == 1) { numComp++; groupComb.push_back(allGroups); }
				
				unweighted = new Unweighted(tmap);
				
			}
			
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "UnifracUnweightedCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void UnifracUnweightedCommand::help(){
	try {
		m->mothurOut("The unifrac.unweighted command can only be executed after a successful read.tree command.\n");
		m->mothurOut("The unifrac.unweighted command parameters are groups, iters, distance, processors and random.  No parameters are required.\n");
		m->mothurOut("The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 1 valid group.\n");
		m->mothurOut("The group names are separated by dashes.  The iters parameter allows you to specify how many random trees you would like compared to your tree.\n");
		m->mothurOut("The distance parameter allows you to create a distance file from the results. The default is false. You may set distance to lt, square or column.\n");
		m->mothurOut("The random parameter allows you to shut off the comparison to random trees. The default is false, meaning compare don't your trees with randomly generated trees.\n");
		m->mothurOut("The processors parameter allows you to specify the number of processors to use. The default is 1.\n");
		m->mothurOut("The unifrac.unweighted command should be in the following format: unifrac.unweighted(groups=yourGroups, iters=yourIters).\n");
		m->mothurOut("Example unifrac.unweighted(groups=A-B-C, iters=500).\n");
		m->mothurOut("The default value for groups is all the groups in your groupfile, and iters is 1000.\n");
		m->mothurOut("The unifrac.unweighted command output two files: .unweighted and .uwsummary their descriptions are in the manual.\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "help");
		exit(1);
	}
}


/***********************************************************/
int UnifracUnweightedCommand::execute() {
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		int start = time(NULL);
		
		userData.resize(numComp,0);  //data[0] = unweightedscore 
		randomData.resize(numComp,0); //data[0] = unweightedscore
		//create new tree with same num nodes and leaves as users
		
		if (numComp < processors) { processors = numComp;  }
		
		outSum << "Tree#" << '\t' << "Groups" << '\t'  <<  "UWScore" <<'\t' << "UWSig" <<  endl;
		m->mothurOut("Tree#\tGroups\tUWScore\tUWSig"); m->mothurOutEndLine();
		
		//get pscores for users trees
		for (int i = 0; i < T.size(); i++) {
			if (m->control_pressed) { 
				outSum.close();
				for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  }
				return 0; 
			}
			
			counter = 0;
			
			if (random)  {  
				output = new ColumnFile(outputDir + m->getSimpleName(globaldata->getTreeFile())  + toString(i+1) + ".unweighted", itersString);
				outputNames.push_back(outputDir + m->getSimpleName(globaldata->getTreeFile())  + toString(i+1) + ".unweighted");
				outputTypes["unweighted"].push_back(outputDir + m->getSimpleName(globaldata->getTreeFile())  + toString(i+1) + ".unweighted");
			}
			
			
			//get unweighted for users tree
			rscoreFreq.resize(numComp);  
			rCumul.resize(numComp);  
			utreeScores.resize(numComp);  
			UWScoreSig.resize(numComp); 

			userData = unweighted->getValues(T[i], processors, outputDir);  //userData[0] = unweightedscore
		
			if (m->control_pressed) { if (random) { delete output;  } outSum.close();  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  }return 0; }
			
			//output scores for each combination
			for(int k = 0; k < numComp; k++) {
				//saves users score
				utreeScores[k].push_back(userData[k]);
				
				//add users score to validscores
				validScores[userData[k]] = userData[k];
			}
		
			//get unweighted scores for random trees - if random is false iters = 0
			for (int j = 0; j < iters; j++) {
		
				//we need a different getValues because when we swap the labels we only want to swap those in each pairwise comparison
				randomData = unweighted->getValues(T[i], "", "", processors, outputDir);
				
				if (m->control_pressed) { if (random) { delete output;  } outSum.close(); for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } return 0; }
			
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
				
				//report progress
//				m->mothurOut("Iter: " + toString(j+1)); m->mothurOutEndLine();	
			}
	
			for(int a = 0; a < numComp; a++) {
				float rcumul = 1.0000;
				
				if (random) {
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
				}else		{	UWScoreSig[a].push_back(0.0);						}
	
			}
			
			if (m->control_pressed) { if (random) { delete output;  } outSum.close(); for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } return 0;  }
			
			//print output files
			printUWSummaryFile(i);
			if (random)  {	printUnweightedFile();	delete output;	}
			if (phylip) {	createPhylipFile(i);		}
			
			rscoreFreq.clear(); 
			rCumul.clear();  
			validScores.clear(); 
			utreeScores.clear();  
			UWScoreSig.clear(); 
		}
		

		outSum.close();
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  }	return 0; }
		
		m->mothurOut("It took " + toString(time(NULL) - start) + " secs to run unifrac.unweighted."); m->mothurOutEndLine();
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "execute");
		exit(1);
	}
}
/***********************************************************/
void UnifracUnweightedCommand::printUnweightedFile() {
	try {
		vector<double> data;
		vector<string> tags;
		
		tags.push_back("Score");
		tags.push_back("RandFreq"); tags.push_back("RandCumul");
			
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
		m->errorOut(e, "UnifracUnweightedCommand", "printUnweightedFile");
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
			m->mothurOut(toString(i+1) + "\t");
			
			if (random) {
				if (UWScoreSig[a][0] > (1/(float)iters)) {
					outSum << setprecision(6) << groupComb[a]  << '\t' << utreeScores[a][0] << '\t' << setprecision(itersString.length()) << UWScoreSig[a][0] << endl;
					cout << setprecision(6)  << groupComb[a]  << '\t' << utreeScores[a][0] << '\t' << setprecision(itersString.length()) << UWScoreSig[a][0] << endl; 
					m->mothurOutJustToLog(groupComb[a]  + "\t" + toString(utreeScores[a][0])  + "\t" + toString(UWScoreSig[a][0])+ "\n"); 
				}else {
					outSum << setprecision(6) << groupComb[a]  << '\t' << utreeScores[a][0] << '\t' << setprecision(itersString.length()) << "<" << (1/float(iters)) << endl;
					cout << setprecision(6)  << groupComb[a]  << '\t' << utreeScores[a][0] << '\t' << setprecision(itersString.length()) << "<" << (1/float(iters)) << endl; 
					m->mothurOutJustToLog(groupComb[a]  + "\t" + toString(utreeScores[a][0])  + "\t<" + toString((1/float(iters))) + "\n"); 
				}
			}else{
				outSum << setprecision(6) << groupComb[a]  << '\t' << utreeScores[a][0] << '\t' << "0.00" << endl;
				cout << setprecision(6)  << groupComb[a]  << '\t' << utreeScores[a][0] << '\t' << "0.00" << endl; 
				m->mothurOutJustToLog(groupComb[a]  + "\t" + toString(utreeScores[a][0])  + "\t0.00\n");
			}
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "printUWSummaryFile");
		exit(1);
	}
}
/***********************************************************/
void UnifracUnweightedCommand::createPhylipFile(int i) {
	try {
		string phylipFileName;
		if ((outputForm == "lt") || (outputForm == "square")) {
			phylipFileName = outputDir + m->getSimpleName(globaldata->getTreeFile())  + toString(i+1) + ".unweighted.phylip.dist";
			outputNames.push_back(phylipFileName); outputTypes["phylip"].push_back(phylipFileName); 
		}else { //column
			phylipFileName = outputDir + m->getSimpleName(globaldata->getTreeFile())  + toString(i+1) + ".unweighted.column.dist";
			outputNames.push_back(phylipFileName); outputTypes["column"].push_back(phylipFileName); 
		}
		
		ofstream out;
		m->openOutputFile(phylipFileName, out);
		
		if ((outputForm == "lt") || (outputForm == "square")) {
			//output numSeqs
			out << globaldata->Groups.size() << endl;
		}
		
		//make matrix with scores in it
		vector< vector<float> > dists;	dists.resize(globaldata->Groups.size());
		for (int i = 0; i < globaldata->Groups.size(); i++) {
			dists[i].resize(globaldata->Groups.size(), 0.0);
		}
		
		//flip it so you can print it
		int count = 0;
		for (int r=0; r<globaldata->Groups.size(); r++) { 
			for (int l = 0; l < r; l++) {
				dists[r][l] = utreeScores[count][0];
				dists[l][r] = utreeScores[count][0];
				count++;
			}
		}
		
		//output to file
		for (int r=0; r<globaldata->Groups.size(); r++) { 
			//output name
			string name = globaldata->Groups[r];
			if (name.length() < 10) { //pad with spaces to make compatible
				while (name.length() < 10) {  name += " ";  }
			}
			
			if (outputForm == "lt") {
				out << name << '\t';
			
				//output distances
				for (int l = 0; l < r; l++) {	out  << dists[r][l] << '\t';  }
				out << endl;
			}else if (outputForm == "square") {
				out << name << '\t';
				
				//output distances
				for (int l = 0; l < globaldata->Groups.size(); l++) {	out << setprecision(11) << dists[r][l] << '\t';  }
				out << endl;
			}else{
				//output distances
				for (int l = 0; l < r; l++) {	
					string otherName = globaldata->Groups[l];
					if (otherName.length() < 10) { //pad with spaces to make compatible
						while (otherName.length() < 10) {  otherName += " ";  }
					}
					
					out  << name << '\t' << otherName << dists[r][l] << endl;  
				}
			}
		}
		out.close();
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "createPhylipFile");
		exit(1);
	}
}
/***********************************************************/



