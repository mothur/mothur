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
vector<string> UnifracUnweightedCommand::setParameters(){	
	try {
		CommandParameter ptree("tree", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(ptree);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pgroup);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pname);
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
		CommandParameter piters("iters", "Number", "", "1000", "", "", "",false,false); parameters.push_back(piters);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "",false,false); parameters.push_back(pprocessors);
		CommandParameter prandom("random", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(prandom);
		CommandParameter pdistance("distance", "Multiple", "column-lt-square", "column", "", "", "",false,false); parameters.push_back(pdistance);
		CommandParameter proot("root", "Boolean", "F", "", "", "", "",false,false); parameters.push_back(proot);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string UnifracUnweightedCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The unifrac.unweighted command parameters are tree, group, name, groups, iters, distance, processors, root and random.  tree parameter is required unless you have valid current tree file.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like analyzed.  You must enter at least 1 valid group.\n";
		helpString += "The group names are separated by dashes.  The iters parameter allows you to specify how many random trees you would like compared to your tree.\n";
		helpString += "The distance parameter allows you to create a distance file from the results. The default is false. You may set distance to lt, square or column.\n";
		helpString += "The random parameter allows you to shut off the comparison to random trees. The default is false, meaning compare don't your trees with randomly generated trees.\n";
		helpString += "The root parameter allows you to include the entire root in your calculations. The default is false, meaning stop at the root for this comparision instead of the root of the entire tree.\n";
		helpString += "The processors parameter allows you to specify the number of processors to use. The default is 1.\n";
		helpString += "The unifrac.unweighted command should be in the following format: unifrac.unweighted(groups=yourGroups, iters=yourIters).\n";
		helpString += "Example unifrac.unweighted(groups=A-B-C, iters=500).\n";
		helpString += "The default value for groups is all the groups in your groupfile, and iters is 1000.\n";
		helpString += "The unifrac.unweighted command output two files: .unweighted and .uwsummary their descriptions are in the manual.\n";
		helpString += "Note: No spaces between parameter labels (i.e. groups), '=' and parameters (i.e.yourGroups).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
UnifracUnweightedCommand::UnifracUnweightedCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
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
/***********************************************************/
UnifracUnweightedCommand::UnifracUnweightedCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
			
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"groups","iters","distance","random","root", "processors","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			map<string,string>::iterator it;
			
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
			
			//check for required parameters
			treefile = validParameter.validFile(parameters, "tree", true);
			if (treefile == "not open") { abort = true; }
			else if (treefile == "not found") { 				//if there is a current design file, use it
				treefile = m->getTreeFile(); 
				if (treefile != "") { m->mothurOut("Using " + treefile + " as input file for the tree parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current tree file and the tree parameter is required."); m->mothurOutEndLine(); abort = true; }								
			}	
			
			//check for required parameters
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }
			else if (groupfile == "not found") { groupfile = ""; }
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }
			else if (namefile == "not found") { namefile = ""; }
			
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				m->splitAtDash(groups, Groups);
				m->Groups = Groups;
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
			
			temp = validParameter.validFile(parameters, "root", false);					if (temp == "not found") { temp = "F"; }
			includeRoot = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			convert(temp, processors); 
			
			if (!random) {  iters = 0;  } //turn off random calcs
			
			//if user selects distance = true and no groups it won't calc the pairwise
			if ((phylip) && (Groups.size() == 0)) {
				groups = "all";
				m->splitAtDash(groups, Groups);
				m->Groups = Groups;
			}
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "UnifracUnweightedCommand", "UnifracUnweightedCommand");
		exit(1);
	}
}

/***********************************************************/
int UnifracUnweightedCommand::execute() {
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
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
					for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } outputTypes.clear();
					m->Groups.clear();
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
	
		sumFile = outputDir + m->getSimpleName(treefile) + ".uwsummary";
		outputNames.push_back(sumFile); outputTypes["uwsummary"].push_back(sumFile);
		m->openOutputFile(sumFile, outSum);
		
		util = new SharedUtil();
		util->setGroups(m->Groups, tmap->namesOfGroups, allGroups, numGroups, "unweighted");	//sets the groups the user wants to analyze
		util->getCombos(groupComb, m->Groups, numComp);
		delete util;
	
		if (numGroups == 1) { numComp++; groupComb.push_back(allGroups); }
		
		unweighted = new Unweighted(tmap, includeRoot);
		
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
				delete tmap; delete unweighted;
				for (int i = 0; i < T.size(); i++) { delete T[i]; }
				outSum.close();
				for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  }
				return 0; 
			}
			
			counter = 0;
			
			if (random)  {  
				output = new ColumnFile(outputDir + m->getSimpleName(treefile)  + toString(i+1) + ".unweighted", itersString);
				outputNames.push_back(outputDir + m->getSimpleName(treefile)  + toString(i+1) + ".unweighted");
				outputTypes["unweighted"].push_back(outputDir + m->getSimpleName(treefile)  + toString(i+1) + ".unweighted");
			}
			
			
			//get unweighted for users tree
			rscoreFreq.resize(numComp);  
			rCumul.resize(numComp);  
			utreeScores.resize(numComp);  
			UWScoreSig.resize(numComp); 

			userData = unweighted->getValues(T[i], processors, outputDir);  //userData[0] = unweightedscore
		
			if (m->control_pressed) { delete tmap; delete unweighted;
				for (int i = 0; i < T.size(); i++) { delete T[i]; }if (random) { delete output;  } outSum.close();  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  }return 0; }
			
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
				
				if (m->control_pressed) { delete tmap; delete unweighted;
					for (int i = 0; i < T.size(); i++) { delete T[i]; }if (random) { delete output;  } outSum.close(); for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } return 0; }
			
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
			
			if (m->control_pressed) { delete tmap; delete unweighted;
				for (int i = 0; i < T.size(); i++) { delete T[i]; }if (random) { delete output;  } outSum.close(); for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } return 0;  }
			
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
		m->Groups.clear();
		delete tmap; delete unweighted;
		for (int i = 0; i < T.size(); i++) { delete T[i]; }
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  }	return 0; }
		
		m->mothurOut("It took " + toString(time(NULL) - start) + " secs to run unifrac.unweighted."); m->mothurOutEndLine();
		
		//set phylip file as new current phylipfile
		string current = "";
		itTypes = outputTypes.find("phylip");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setPhylipFile(current); }
		}
		
		//set column file as new current columnfile
		itTypes = outputTypes.find("column");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setColumnFile(current); }
		}
		
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
			phylipFileName = outputDir + m->getSimpleName(treefile)  + toString(i+1) + ".unweighted.phylip.dist";
			outputNames.push_back(phylipFileName); outputTypes["phylip"].push_back(phylipFileName); 
		}else { //column
			phylipFileName = outputDir + m->getSimpleName(treefile)  + toString(i+1) + ".unweighted.column.dist";
			outputNames.push_back(phylipFileName); outputTypes["column"].push_back(phylipFileName); 
		}
		
		ofstream out;
		m->openOutputFile(phylipFileName, out);
		
		if ((outputForm == "lt") || (outputForm == "square")) {
			//output numSeqs
			out << m->Groups.size() << endl;
		}
		
		//make matrix with scores in it
		vector< vector<float> > dists;	dists.resize(m->Groups.size());
		for (int i = 0; i < m->Groups.size(); i++) {
			dists[i].resize(m->Groups.size(), 0.0);
		}
		
		//flip it so you can print it
		int count = 0;
		for (int r=0; r<m->Groups.size(); r++) { 
			for (int l = 0; l < r; l++) {
				dists[r][l] = utreeScores[count][0];
				dists[l][r] = utreeScores[count][0];
				count++;
			}
		}
		
		//output to file
		for (int r=0; r<m->Groups.size(); r++) { 
			//output name
			string name = m->Groups[r];
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
				for (int l = 0; l < m->Groups.size(); l++) {	out << dists[r][l] << '\t';  }
				out << endl;
			}else{
				//output distances
				for (int l = 0; l < r; l++) {	
					string otherName = m->Groups[l];
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
}/*****************************************************************/
int UnifracUnweightedCommand::readNamesFile() {
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
		m->errorOut(e, "UnifracUnweightedCommand", "readNamesFile");
		exit(1);
	}
}
/***********************************************************/




