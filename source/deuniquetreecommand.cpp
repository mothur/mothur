/*
 *  deuniquetreecommand.cpp
 *  Mothur
 *
 *  Created by westcott on 5/27/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "deuniquetreecommand.h"
#include "treereader.h"

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
		
		TreeReader* reader = new TreeReader(treefile, "", namefile);
        vector<Tree*> T = reader->getTrees();
        map<string, string> nameMap = reader->getNameMap();
        delete reader;		
		
		//print new Tree
		string outputFile = outputDir + m->getRootName(m->getSimpleName(treefile)) + "deunique.tre";
		outputNames.push_back(outputFile); outputTypes["tree"].push_back(outputFile);
		ofstream out;
		m->openOutputFile(outputFile, out);
		T[0]->print(out, nameMap);
		out.close();
		
        delete (T[0]->getTreeMap());
		for (int i = 0; i < T.size(); i++) { delete T[i]; }
				
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
/***********************************************************/




