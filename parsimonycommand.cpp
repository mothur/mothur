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
vector<string> ParsimonyCommand::setParameters(){	
	try {
		CommandParameter ptree("tree", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(ptree);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pgroup);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pname);
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
		CommandParameter prandom("random", "String", "", "", "", "", "",false,false); parameters.push_back(prandom);
		CommandParameter piters("iters", "Number", "", "1000", "", "", "",false,false); parameters.push_back(piters);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "",false,false); parameters.push_back(pprocessors);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ParsimonyCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The parsimony command parameters are tree, group, name, random, groups, processors and iters.  tree parameter is required unless you have valid current tree file or are using random.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 1 valid group.\n";
		helpString += "The group names are separated by dashes.  The iters parameter allows you to specify how many random trees you would like compared to your tree.\n";
		helpString += "The parsimony command should be in the following format: parsimony(random=yourOutputFilename, groups=yourGroups, iters=yourIters).\n";
		helpString += "The processors parameter allows you to specify the number of processors to use. The default is 1.\n";
		helpString += "Example parsimony(random=out, iters=500).\n";
		helpString += "The default value for random is "" (meaning you want to use the trees in your inputfile, randomtree=out means you just want the random distribution of trees outputted to out.rd_parsimony),\n";
		helpString += "and iters is 1000.  The parsimony command output two files: .parsimony and .psummary their descriptions are in the manual.\n";
		helpString += "Note: No spaces between parameter labels (i.e. random), '=' and parameters (i.e.yourOutputFilename).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
ParsimonyCommand::ParsimonyCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["parsimony"] = tempOutNames;
		outputTypes["psummary"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "ParsimonyCommand");
		exit(1);
	}
}
/***********************************************************/
ParsimonyCommand::ParsimonyCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		Groups.clear();
			
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			map<string,string>::iterator it;
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["parsimony"] = tempOutNames;
			outputTypes["psummary"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("tree");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["tree"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
			}
			
			m->runParse = true;
			m->clearGroups();
			m->clearAllGroups();
			m->Treenames.clear();
			m->names.clear();
			
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}
			
			randomtree = validParameter.validFile(parameters, "random", false);		if (randomtree == "not found") { randomtree = ""; }
			
			//are you trying to use parsimony without reading a tree or saying you want random distribution
			if (randomtree == "")  {
				//check for required parameters
				treefile = validParameter.validFile(parameters, "tree", true);
				if (treefile == "not open") { treefile = ""; abort = true; }
				else if (treefile == "not found") { 				//if there is a current design file, use it
					treefile = m->getTreeFile(); 
					if (treefile != "") { m->mothurOut("Using " + treefile + " as input file for the tree parameter."); m->mothurOutEndLine(); }
					else { 	m->mothurOut("You have no current tree file and the tree parameter is required."); m->mothurOutEndLine(); abort = true; }								
				}else { m->setTreeFile(treefile); }	
				
				//check for required parameters
				groupfile = validParameter.validFile(parameters, "group", true);
				if (groupfile == "not open") { abort = true; }
				else if (groupfile == "not found") { groupfile = ""; }
				else { m->setGroupFile(groupfile); }
				
				namefile = validParameter.validFile(parameters, "name", true);
				if (namefile == "not open") { namefile = ""; abort = true; }
				else if (namefile == "not found") { namefile = ""; }
				else { m->setNameFile(namefile); }
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			string outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	if (randomtree == "")  { outputDir += m->hasPath(treefile); } }
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; m->clearGroups(); }
			else { 
				m->splitAtDash(groups, Groups);
				m->setGroups(Groups);
			}
				
			itersString = validParameter.validFile(parameters, "iters", false);			if (itersString == "not found") { itersString = "1000"; }
			m->mothurConvert(itersString, iters); 
			
			string temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors);
			
			if (namefile == "") {
				vector<string> files; files.push_back(treefile);
				parser.getNameFile(files);
			}
			
		}

	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "ParsimonyCommand");
		exit(1);
	}
}
/***********************************************************/
int ParsimonyCommand::execute() {
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		
		//randomtree will tell us if user had their own treefile or if they just want the random distribution
		//user has entered their own tree
		if (randomtree == "") { 
			
			m->setTreeFile(treefile);
			
			if (groupfile != "") {
				//read in group map info.
				tmap = new TreeMap(groupfile);
				tmap->readMap();
			}else{ //fake out by putting everyone in one group
				Tree* tree = new Tree(treefile); delete tree;  //extracts names from tree to make faked out groupmap
				tmap = new TreeMap();
				
				for (int i = 0; i < m->Treenames.size(); i++) { tmap->addSeq(m->Treenames[i], "Group1"); }
			}
			
			if (namefile != "") { readNamesFile(); }
			
			read = new ReadNewickTree(treefile);
			int readOk = read->read(tmap); 
			
			if (readOk != 0) { m->mothurOut("Read Terminated."); m->mothurOutEndLine(); delete tmap; delete read; return 0; }
			
			read->AssembleTrees();
			T = read->getTrees();
			delete read;

			//make sure all files match
			//if you provide a namefile we will use the numNames in the namefile as long as the number of unique match the tree names size.
			int numNamesInTree;
			if (namefile != "")  {  
				if (numUniquesInName == m->Treenames.size()) {  numNamesInTree = nameMap.size();  }
				else {   numNamesInTree = m->Treenames.size();  }
			}else {  numNamesInTree = m->Treenames.size();  }
			
			
			//output any names that are in group file but not in tree
			if (numNamesInTree < tmap->getNumSeqs()) {
				for (int i = 0; i < tmap->namesOfSeqs.size(); i++) {
					//is that name in the tree?
					int count = 0;
					for (int j = 0; j < m->Treenames.size(); j++) {
						if (tmap->namesOfSeqs[i] == m->Treenames[j]) { break; } //found it
						count++;
					}
					
					if (m->control_pressed) { 
						delete tmap; for (int i = 0; i < T.size(); i++) { delete T[i]; }
						for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } outputTypes.clear();
						m->clearGroups();
						return 0;
					}
					
					//then you did not find it so report it 
					if (count == m->Treenames.size()) { 
						//if it is in your namefile then don't remove
						map<string, string>::iterator it = nameMap.find(tmap->namesOfSeqs[i]);
						
						if (it == nameMap.end()) {
							m->mothurOut(tmap->namesOfSeqs[i] + " is in your groupfile and not in your tree. It will be disregarded."); m->mothurOutEndLine();
							tmap->removeSeq(tmap->namesOfSeqs[i]);
							i--; //need this because removeSeq removes name from namesOfSeqs
						}
					}
				}
			}
				
			if(outputDir == "") { outputDir += m->hasPath(treefile); }
			output = new ColumnFile(outputDir + m->getSimpleName(treefile)  +  ".parsimony", itersString);
			outputNames.push_back(outputDir + m->getSimpleName(treefile)  +  ".parsimony");
			outputTypes["parsimony"].push_back(outputDir + m->getSimpleName(treefile)  +  ".parsimony");
				
			sumFile = outputDir + m->getSimpleName(treefile) + ".psummary";
			m->openOutputFile(sumFile, outSum);
			outputNames.push_back(sumFile);
			outputTypes["psummary"].push_back(sumFile);
		}else { //user wants random distribution
			getUserInput();
				
			if(outputDir == "") { outputDir += m->hasPath(randomtree); }
			output = new ColumnFile(outputDir+ m->getSimpleName(randomtree), itersString);
			outputNames.push_back(outputDir+ m->getSimpleName(randomtree));
			outputTypes["parsimony"].push_back(outputDir+ m->getSimpleName(randomtree));
		}
			
		//set users groups to analyze
		util = new SharedUtil();
		vector<string> mGroups = m->getGroups();
		vector<string> tGroups = tmap->getNamesOfGroups();
		util->setGroups(mGroups, tGroups, allGroups, numGroups, "parsimony");	//sets the groups the user wants to analyze
		util->getCombos(groupComb, mGroups, numComp);
		m->setGroups(mGroups);
		delete util;
			
		if (numGroups == 1) { numComp++; groupComb.push_back(allGroups); }
			
		pars = new Parsimony(tmap);
		counter = 0;
	
		Progress* reading;
		reading = new Progress("Comparing to random:", iters);
		
		if (m->control_pressed) { 
			delete reading; delete pars; delete output;
			delete tmap; for (int i = 0; i < T.size(); i++) { delete T[i]; }
			if (randomtree == "") {  outSum.close();  }
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } outputTypes.clear();
			m->clearGroups();
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
					delete reading; delete pars; delete output;
					delete tmap; for (int i = 0; i < T.size(); i++) { delete T[i]; }
					if (randomtree == "") {  outSum.close();  }
					for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } outputTypes.clear();
					m->clearGroups();
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
				randT = new Tree(tmap);

				//create random relationships between nodes
				randT->assembleRandomTree();

				//get pscore of random tree
				randomData = pars->getValues(randT, processors, outputDir);
				
				if (m->control_pressed) { 
					delete reading; delete pars; delete output; delete randT;
					if (randomtree == "") {  outSum.close();  }
					for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } outputTypes.clear();
					delete tmap; for (int i = 0; i < T.size(); i++) { delete T[i]; }
					m->clearGroups();
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
				randT = new Tree(tmap);
				//create random relationships between nodes

				randT->assembleRandomTree();
				
				if (m->control_pressed) { 
					delete reading; delete pars; delete output; delete randT;
					delete tmap; 
					for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } outputTypes.clear();
					m->clearGroups();
					return 0;
				}


				//get pscore of random tree
				randomData = pars->getValues(randT, processors, outputDir);
				
				if (m->control_pressed) { 
					delete reading; delete pars;  delete output; delete randT;
					delete tmap; 
					for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } outputTypes.clear();
					m->clearGroups();
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
				delete reading; delete pars; delete output;
				delete tmap; for (int i = 0; i < T.size(); i++) { delete T[i]; }
				if (randomtree == "") {  outSum.close();  }
				for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } outputTypes.clear();
				m->clearGroups();
				return 0;
		}
		
		//finish progress bar
		reading->finish();
		delete reading;

		
		printParsimonyFile();
		if (randomtree == "") { printUSummaryFile(); }
		
		//reset groups parameter
		m->clearGroups(); 
		
		delete pars; delete output; 
		delete tmap; for (int i = 0; i < T.size(); i++) { delete T[i]; }
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } outputTypes.clear(); return 0;}
		
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
			tmap->addGroup(toString(i));
			
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
		
		m->Treenames = tmap->namesOfSeqs; 
		
	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "getUserInput");
		exit(1);
	}
}
/*****************************************************************/
int ParsimonyCommand::readNamesFile() {
	try {
		m->names.clear();
		numUniquesInName = 0;
		
		ifstream in;
		m->openInputFile(namefile, in);
		
		string first, second;
		map<string, string>::iterator itNames;
		
		while(!in.eof()) {
			in >> first >> second; m->gobble(in);
			
			numUniquesInName++;
			
			itNames = m->names.find(first);
			if (itNames == m->names.end()) {  
				m->names[first] = second; 
				
				//we need a list of names in your namefile to use above when removing extra seqs above so we don't remove them
				vector<string> dupNames;
				m->splitAtComma(second, dupNames);
				
				for (int i = 0; i < dupNames.size(); i++) {	
					nameMap[dupNames[i]] = dupNames[i]; 
					if ((groupfile == "") && (i != 0)) { tmap->addSeq(dupNames[i], "Group1"); } 
				}
			}else {  m->mothurOut(first + " has already been seen in namefile, disregarding names file."); m->mothurOutEndLine(); in.close(); m->names.clear(); namefile = ""; return 1; }			
		}
		in.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ParsimonyCommand", "readNamesFile");
		exit(1);
	}
}
/***********************************************************/


