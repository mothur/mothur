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
        
        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["tree"] = tempOutNames;
		
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
/***********************************************************/
DeuniqueTreeCommand::DeuniqueTreeCommand(string option)  {
	try {
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			treefile = validParameter.validFile(parameters, "tree");
			if (treefile == "not open") { abort = true; }
			else if (treefile == "not found") { 				//if there is a current design file, use it
				treefile = current->getTreeFile(); 
				if (treefile != "") { m->mothurOut("Using " + treefile + " as input file for the tree parameter.\n");  }
				else { 	m->mothurOut("You have no current tree file and the tree parameter is required.\n");  abort = true; }
			}else { current->setTreeFile(treefile); }	
			
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { abort = true; }
			else if (namefile == "not found") { 				//if there is a current design file, use it
				namefile = current->getNameFile(); 
				if (namefile != "") { m->mothurOut("Using " + namefile + " as input file for the name parameter.\n");  }
				else { 	m->mothurOut("You have no current name file and the name parameter is required.\n");  abort = true; }
			}else { current->setNameFile(namefile); }
			
					if (outputdir == ""){    outputdir = util.hasPath(treefile);	}
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
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(treefile));
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
		
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "DeuniqueTreeCommand", "execute");
		exit(1);
	}
}
/***********************************************************/




