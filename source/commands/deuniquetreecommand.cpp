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
		CommandParameter ptree("tree", "InputTypes", "", "", "none", "none", "none","tree",false,true,true); parameters.push_back(ptree);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(pname);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
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
string DeuniqueTreeCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "tree") {  pattern = "[filename],deunique.tre"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "DeuniqueTreeCommand", "getOutputPattern");
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
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("tree");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["tree"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
			}
			
            //check for required parameters
			treefile = validParameter.validFile(parameters, "tree");
			if (treefile == "not open") { abort = true; }
			else if (treefile == "not found") { 				//if there is a current design file, use it
				treefile = current->getTreeFile(); 
				if (treefile != "") { m->mothurOut("Using " + treefile + " as input file for the tree parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current tree file and the tree parameter is required."); m->mothurOutEndLine(); abort = true; }								
			}else { current->setTreeFile(treefile); }	
			
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { abort = true; }
			else if (namefile == "not found") { 				//if there is a current design file, use it
				namefile = current->getNameFile(); 
				if (namefile != "") { m->mothurOut("Using " + namefile + " as input file for the name parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current name file and the name parameter is required."); m->mothurOutEndLine(); abort = true; }								
			}else { current->setNameFile(namefile); }
			
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = util.hasPath(treefile);	}
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
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		current->setTreeFile(treefile);
		
		TreeReader* reader = new TreeReader(treefile, "", namefile);
        vector<Tree*> T = reader->getTrees();
        map<string, string> nameMap;
        util.readNames(namefile, nameMap);
        delete reader;		
		
		//print new Tree
        map<string, string> variables; 
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(treefile));
		string outputFile = getOutputFileName("tree", variables);
		outputNames.push_back(outputFile); outputTypes["tree"].push_back(outputFile);
		ofstream out;
		util.openOutputFile(outputFile, out);
		T[0]->print(out, nameMap);
		out.close();
		
        delete (T[0]->getCountTable());
		for (int i = 0; i < T.size(); i++) { delete T[i]; }
				
		//set phylip file as new current phylipfile
		string currentName = "";
		itTypes = outputTypes.find("tree");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setTreeFile(currentName); }
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




