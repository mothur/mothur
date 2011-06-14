/*
 *  deuniquetreecommand.cpp
 *  Mothur
 *
 *  Created by westcott on 5/27/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "deuniquetreecommand.h"

//**********************************************************************************************************************
vector<string> DeuniqueTreeCommand::setParameters(){	
	try {
		CommandParameter ptree("tree", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(ptree);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pname);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "DeuniqueTreeCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string DeuniqueTreeCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The deunique.tree command parameters are tree and name.  Both parameters are required unless you have valid current files.\n";
		helpString += "The deunique.tree command should be in the following format: deunique.tree(tree=yourTreeFile, name=yourNameFile).\n";
		helpString += "Example deunique.tree(tree=abrecovery.tree, name=abrecovery.names).\n";
		helpString += "Note: No spaces between parameter labels (i.e. tree), '=' and parameters (i.e.yourTreeFile).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "DeuniqueTreeCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
DeuniqueTreeCommand::DeuniqueTreeCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["tree"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "DeuniqueTreeCommand", "DeuniqueTreeCommand");
		exit(1);
	}
}
/***********************************************************/
DeuniqueTreeCommand::DeuniqueTreeCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
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
			outputTypes["tree"] = tempOutNames;
						
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
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
			}
			
			m->runParse = true;
			m->Groups.clear();
			m->namesOfGroups.clear();
			m->Treenames.clear();
			m->names.clear();
			
			//check for required parameters
			treefile = validParameter.validFile(parameters, "tree", true);
			if (treefile == "not open") { abort = true; }
			else if (treefile == "not found") { 				//if there is a current design file, use it
				treefile = m->getTreeFile(); 
				if (treefile != "") { m->mothurOut("Using " + treefile + " as input file for the tree parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current tree file and the tree parameter is required."); m->mothurOutEndLine(); abort = true; }								
			}else { m->setTreeFile(treefile); }	
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }
			else if (namefile == "not found") { 				//if there is a current design file, use it
				namefile = m->getNameFile(); 
				if (namefile != "") { m->mothurOut("Using " + namefile + " as input file for the name parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current name file and the name parameter is required."); m->mothurOutEndLine(); abort = true; }								
			}else { m->setNameFile(namefile); }
			
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = m->hasPath(treefile);	}
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "DeuniqueTreeCommand", "DeuniqueTreeCommand");
		exit(1);
	}
}

/***********************************************************/
int DeuniqueTreeCommand::execute() {
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		m->setTreeFile(treefile);
		
		//extracts names from tree to make faked out groupmap
		Tree* tree = new Tree(treefile); delete tree;  
		tmap = new TreeMap();
		for (int i = 0; i < m->Treenames.size(); i++) { tmap->addSeq(m->Treenames[i], "Group1"); }
		
		if (m->control_pressed) {  delete tmap;  return 0; }
		
		readNamesFile(); 
		
		if (m->control_pressed) {  delete tmap;  return 0; }
		
		ReadTree* read = new ReadNewickTree(treefile);
		int readOk = read->read(tmap); 
		if (readOk != 0) { m->mothurOut("Read Terminated."); m->mothurOutEndLine(); delete tmap; delete read; return 0; }
		
		read->AssembleTrees();
		vector<Tree*> T = read->getTrees();
		delete read;
		
		//make sure all files match
		//if you provide a namefile we will use the numNames in the namefile as long as the number of unique match the tree names size.
		int numNamesInTree;
		if (numUniquesInName == m->Treenames.size()) {  numNamesInTree = nameMap.size();  }
		else {   numNamesInTree = m->Treenames.size();  }
		
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
		
		
		//print new Tree
		string outputFile = outputDir + m->getRootName(m->getSimpleName(treefile)) + "deunique.tre";
		outputNames.push_back(outputFile); outputTypes["tree"].push_back(outputFile);
		ofstream out;
		m->openOutputFile(outputFile, out);
		T[0]->print(out, "deunique");
		out.close();
		
		delete tmap; for (int i = 0; i < T.size(); i++) { delete T[i]; }
				
		//set phylip file as new current phylipfile
		string current = "";
		itTypes = outputTypes.find("tree");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setTreeFile(current); }
		}
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "DeuniqueTreeCommand", "execute");
		exit(1);
	}
}
/*****************************************************************/
int DeuniqueTreeCommand::readNamesFile() {
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
					if (i != 0) { tmap->addSeq(dupNames[i], "Group1"); } 
				}
			}else {  m->mothurOut(first + " has already been seen in namefile, aborting."); m->mothurOutEndLine(); in.close(); m->names.clear(); m->control_pressed = true; return 1; }			
		}
		in.close();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "DeuniqueTreeCommand", "readNamesFile");
		exit(1);
	}
}
/***********************************************************/




