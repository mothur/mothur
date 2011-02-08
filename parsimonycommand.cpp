/*
 *  parsimonycommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/26/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "parsimonycommand.h"

//**********************************************************************************************************************
vector<string> ParsimonyCommand::getValidParameters(){	
	try {
		string Array[] =  {"random","groups","iters","processors","outputdir","inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
ParsimonyCommand::ParsimonyCommand(){	
	try {
		abort = true; calledHelp = true; 
		vector<string> tempOutNames;
		outputTypes["parsimony"] = tempOutNames;
		outputTypes["psummary"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "ParsimonyCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> ParsimonyCommand::getRequiredParameters(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> ParsimonyCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "getRequiredFiles");
		exit(1);
	}
}
/***********************************************************/
ParsimonyCommand::ParsimonyCommand(string option)  {
	try {
		globaldata = GlobalData::getInstance();
		abort = false; calledHelp = false;   
		Groups.clear();
			
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"random","groups","processors","iters","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["parsimony"] = tempOutNames;
			outputTypes["psummary"] = tempOutNames;
			
			randomtree = validParameter.validFile(parameters, "random", false);		if (randomtree == "not found") { randomtree = ""; }
			
			//are you trying to use parsimony without reading a tree or saying you want random distribution
			if (randomtree == "")  {
				if (globaldata->gTree.size() == 0) {
					m->mothurOut("You must read a treefile and a groupfile or set the randomtree parameter to the output filename you wish, before you may execute the parsimony command."); m->mothurOutEndLine(); abort = true;  }
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			string outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	if (randomtree == "")  { outputDir += m->hasPath(globaldata->inputFileName); } }
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; globaldata->Groups.clear(); }
			else { 
				m->splitAtDash(groups, Groups);
				globaldata->Groups = Groups;
			}
				
			itersString = validParameter.validFile(parameters, "iters", false);			if (itersString == "not found") { itersString = "1000"; }
			convert(itersString, iters); 
			
			string temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = "1";				}
			convert(temp, processors); 
						
			if (abort == false) {
				//randomtree will tell us if user had their own treefile or if they just want the random distribution
				//user has entered their own tree
				if (randomtree == "") { 
					T = globaldata->gTree;
					tmap = globaldata->gTreemap;
					
					if(outputDir == "") { outputDir += m->hasPath(globaldata->getTreeFile()); }
					output = new ColumnFile(outputDir + m->getSimpleName(globaldata->getTreeFile())  +  ".parsimony", itersString);
					outputNames.push_back(outputDir + m->getSimpleName(globaldata->getTreeFile())  +  ".parsimony");
					outputTypes["parsimony"].push_back(outputDir + m->getSimpleName(globaldata->getTreeFile())  +  ".parsimony");
					
					sumFile = outputDir + m->getSimpleName(globaldata->getTreeFile()) + ".psummary";
					m->openOutputFile(sumFile, outSum);
					outputNames.push_back(sumFile);
					outputTypes["psummary"].push_back(sumFile);
				}else { //user wants random distribution
					savetmap = globaldata->gTreemap;
					getUserInput();
					
					if(outputDir == "") { outputDir += m->hasPath(randomtree); }
					output = new ColumnFile(outputDir+ m->getSimpleName(randomtree), itersString);
					outputNames.push_back(outputDir+ m->getSimpleName(randomtree));
					outputTypes["parsimony"].push_back(outputDir+ m->getSimpleName(randomtree));
				}
				
				//set users groups to analyze
				util = new SharedUtil();
				util->setGroups(globaldata->Groups, tmap->namesOfGroups, allGroups, numGroups, "parsimony");	//sets the groups the user wants to analyze
				util->getCombos(groupComb, globaldata->Groups, numComp);
				
				if (numGroups == 1) { numComp++; groupComb.push_back(allGroups); }
				
				pars = new Parsimony(tmap);
				counter = 0;
				
			}
			
		}

	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "ParsimonyCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void ParsimonyCommand::help(){
	try {
		m->mothurOut("The parsimony command can only be executed after a successful read.tree command, unless you use the random parameter.\n");
		m->mothurOut("The parsimony command parameters are random, groups, processors and iters.  No parameters are required.\n");
		m->mothurOut("The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 1 valid group.\n");
		m->mothurOut("The group names are separated by dashes.  The iters parameter allows you to specify how many random trees you would like compared to your tree.\n");
		m->mothurOut("The parsimony command should be in the following format: parsimony(random=yourOutputFilename, groups=yourGroups, iters=yourIters).\n");
		m->mothurOut("The processors parameter allows you to specify the number of processors to use. The default is 1.\n");
		m->mothurOut("Example parsimony(random=out, iters=500).\n");
		m->mothurOut("The default value for random is "" (meaning you want to use the trees in your inputfile, randomtree=out means you just want the random distribution of trees outputted to out.rd_parsimony),\n");
		m->mothurOut("and iters is 1000.  The parsimony command output two files: .parsimony and .psummary their descriptions are in the manual.\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. random), '=' and parameters (i.e.yourOutputFilename).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "help");
		exit(1);
	}
}


/***********************************************************/
int ParsimonyCommand::execute() {
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
	
		Progress* reading;
		reading = new Progress("Comparing to random:", iters);
		
		if (m->control_pressed) { 
			delete reading; delete pars; delete util; delete output;
			if (randomtree == "") {  outSum.close();  }
			for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } outputTypes.clear();
			globaldata->Groups.clear();
			return 0;
		}
			
		
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
				userData = pars->getValues(T[i], processors, outputDir);  //data = AB, AC, BC, ABC.
				
				if (m->control_pressed) { 
					delete reading; delete pars; delete util; delete output;
					if (randomtree == "") {  outSum.close();  }
					for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } outputTypes.clear();
					globaldata->Groups.clear();
					return 0;
				}


				//output scores for each combination
				for(int k = 0; k < numComp; k++) {

					//update uscoreFreq
					map<int,double>::iterator it = uscoreFreq[k].find(userData[k]);
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
				randomData = pars->getValues(randT, processors, outputDir);
				
				if (m->control_pressed) { 
					delete reading; delete pars; delete util; delete output; delete randT;
					if (randomtree == "") {  outSum.close();  }
					for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } outputTypes.clear();
					globaldata->Groups.clear();
					return 0;
				}
					
				for(int r = 0; r < numComp; r++) {
					//add trees pscore to map of scores
					map<int,double>::iterator it = rscoreFreq[r].find(randomData[r]);
					if (it != rscoreFreq[r].end()) {//already have that score
						rscoreFreq[r][randomData[r]]++;
					}else{//first time we have seen this score
						rscoreFreq[r][randomData[r]] = 1;
					}
			
					//add randoms score to validscores
					validScores[randomData[r]] = randomData[r];
				}
				
				//update progress bar
				reading->update(j);
				
				delete randT;
			}

		}else {
			//get pscores for random trees
			for (int j = 0; j < iters; j++) {
								
				//create new tree with same num nodes and leaves as users
				randT = new Tree();
				//create random relationships between nodes

				randT->assembleRandomTree();
				
				if (m->control_pressed) { 
					delete reading; delete pars; delete util; delete output; delete randT;
					globaldata->gTreemap = savetmap; 
					for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } outputTypes.clear();
					globaldata->Groups.clear();
					return 0;
				}


				//get pscore of random tree
				randomData = pars->getValues(randT, processors, outputDir);
				
				if (m->control_pressed) { 
					delete reading; delete pars; delete util; delete output; delete randT;
					globaldata->gTreemap = savetmap; 
					for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } outputTypes.clear();
					globaldata->Groups.clear();
					return 0;
				}
			
				for(int r = 0; r < numComp; r++) {
					//add trees pscore to map of scores
					map<int,double>::iterator it = rscoreFreq[r].find(randomData[r]);
					if (it != rscoreFreq[r].end()) {//already have that score
						rscoreFreq[r][randomData[r]]++;
					}else{//first time we have seen this score
						rscoreFreq[r][randomData[r]] = 1;
					}
			
					//add randoms score to validscores
					validScores[randomData[r]] = randomData[r];
				}
				
				//update progress bar
				reading->update(j);
				
				delete randT;
			}
		}

		for(int a = 0; a < numComp; a++) {
			float rcumul = 0.0000;
			float ucumul = 0.0000;
			//this loop fills the cumulative maps and put 0.0000 in the score freq map to make it easier to print.
			for (map<int,double>::iterator it = validScores.begin(); it != validScores.end(); it++) { 
				if (randomtree == "") {
					map<int,double>::iterator it2 = uscoreFreq[a].find(it->first);
					//user data has that score 
					if (it2 != uscoreFreq[a].end()) { uscoreFreq[a][it->first] /= T.size(); ucumul+= it2->second;  }
					else { uscoreFreq[a][it->first] = 0.0000; } //no user trees with that score
					//make uCumul map
					uCumul[a][it->first] = ucumul;
				}
			
				//make rscoreFreq map and rCumul
				map<int,double>::iterator it2 = rscoreFreq[a].find(it->first);
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
		
		if (m->control_pressed) { 
				delete reading; delete pars; delete util; delete output;
				if (randomtree == "") {  outSum.close();  }
				else { globaldata->gTreemap = savetmap; }
				for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } outputTypes.clear();
				globaldata->Groups.clear();
				return 0;
		}
		
		//finish progress bar
		reading->finish();
		delete reading;

		
		printParsimonyFile();
		if (randomtree == "") { printUSummaryFile(); }
		
		//reset globaldata's treemap if you just did random distrib
		if (randomtree != "") {
			//memory leak prevention
			//if (globaldata->gTreemap != NULL) { delete globaldata->gTreemap;  }
			globaldata->gTreemap = savetmap;
		}
		
		//reset groups parameter
		globaldata->Groups.clear(); 
		
		if (m->control_pressed) { 
			delete pars; delete util; delete output;
			for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } outputTypes.clear();
			return 0;
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "execute");
		exit(1);
	}
}

/***********************************************************/
void ParsimonyCommand::printParsimonyFile() {
	try {
		vector<double> data;
		vector<string> tags;
		
		if (randomtree == "") {
			tags.push_back("Score"); tags.push_back("UserFreq"); tags.push_back("UserCumul"); tags.push_back("RandFreq"); tags.push_back("RandCumul");
		}else {
			tags.push_back("Score"); tags.push_back("RandFreq"); tags.push_back("RandCumul");
		}

		for(int a = 0; a < numComp; a++) {
			output->initFile(groupComb[a], tags);
			//print each line
			for (map<int,double>::iterator it = validScores.begin(); it != validScores.end(); it++) { 
				if (randomtree == "") {
					data.push_back(it->first);  data.push_back(uscoreFreq[a][it->first]); data.push_back(uCumul[a][it->first]); data.push_back(rscoreFreq[a][it->first]); data.push_back(rCumul[a][it->first]); 
				}else{
					data.push_back(it->first);  data.push_back(rscoreFreq[a][it->first]); data.push_back(rCumul[a][it->first]); 
				}
				output->output(data);
				data.clear();
			} 
			output->resetFile();
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "printParsimonyFile");
		exit(1);
	}
}
/***********************************************************/
int ParsimonyCommand::printUSummaryFile() {
	try {
		//column headers
		outSum << "Tree#" << '\t' << "Groups" << '\t'  <<  "ParsScore" << '\t' << "ParsSig" <<  endl;
		m->mothurOut("Tree#\tGroups\tParsScore\tParsSig"); m->mothurOutEndLine();
		
		//format output
		outSum.setf(ios::fixed, ios::floatfield); outSum.setf(ios::showpoint);
		
		
		//print each line
		for (int i = 0; i< T.size(); i++) {
			for(int a = 0; a < numComp; a++) {
				if (m->control_pressed) {  outSum.close(); return 0; }
				if (UScoreSig[a][i] > (1/(float)iters)) {
					outSum << setprecision(6) << i+1 << '\t' << groupComb[a]  << '\t' << userTreeScores[a][i] << setprecision(itersString.length()) << '\t' << UScoreSig[a][i] << endl;
					cout << setprecision(6) << i+1 << '\t' << groupComb[a]  << '\t' << userTreeScores[a][i] << setprecision(itersString.length()) << '\t' << UScoreSig[a][i] << endl;
					m->mothurOutJustToLog(toString(i+1) + "\t" + groupComb[a] + "\t" + toString(userTreeScores[a][i]) + "\t" + toString(UScoreSig[a][i])); m->mothurOutEndLine();
				}else {
					outSum << setprecision(6) << i+1 << '\t' << groupComb[a] << '\t' << userTreeScores[a][i] << setprecision(itersString.length())  << '\t' << "<" << (1/float(iters)) << endl;
					cout << setprecision(6) << i+1 << '\t' << groupComb[a] << '\t' << userTreeScores[a][i] << setprecision(itersString.length()) << '\t' << "<" << (1/float(iters)) << endl;
					m->mothurOutJustToLog(toString(i+1) + "\t" + groupComb[a] + "\t" + toString(userTreeScores[a][i]) + "\t" + toString((1/float(iters)))); m->mothurOutEndLine();
				}
			}
		}
		
		outSum.close();
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "printUSummaryFile");
		exit(1);
	}
}

/***********************************************************/
void ParsimonyCommand::getUserInput() {
	try {
	
		//create treemap
		tmap = new TreeMap();

		m->mothurOut("Please enter the number of groups you would like to analyze: ");
		cin >> numGroups;
		m->mothurOutJustToLog(toString(numGroups)); m->mothurOutEndLine();
				
		int num, count;
		count = 1;
		numEachGroup.resize(numGroups, 0);  
		
		for (int i = 1; i <= numGroups; i++) {
			m->mothurOut("Please enter the number of sequences in group " + toString(i) +  ": ");
			cin >> num;
			m->mothurOutJustToLog(toString(num)); m->mothurOutEndLine();
				
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
		//if (globaldata->gTreemap != NULL) { delete globaldata->gTreemap;  }
		globaldata->gTreemap = tmap;
		globaldata->Treenames = tmap->namesOfSeqs; 
		
	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "getUserInput");
		exit(1);
	}
}

/***********************************************************/


