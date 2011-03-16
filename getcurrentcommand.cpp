/*
 *  getcurrentcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 3/16/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "getcurrentcommand.h"


//**********************************************************************************************************************
vector<string> GetCurrentCommand::getValidParameters(){	
	try {
		string Array[] =  {"outputdir","inputdir","clear"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetCurrentCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
GetCurrentCommand::GetCurrentCommand(){	
	try {
		abort = true; calledHelp = true; 
	}
	catch(exception& e) {
		m->errorOut(e, "GetCurrentCommand", "GetCurrentCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> GetCurrentCommand::getRequiredParameters(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetCurrentCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> GetCurrentCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetCurrentCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************
GetCurrentCommand::GetCurrentCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"outputdir","inputdir","clear"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			clearTypes = validParameter.validFile(parameters, "clear", false);			
			if (clearTypes == "not found") { clearTypes = ""; }
			else { m->splitAtDash(clearTypes, types);	}
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "GetCurrentCommand", "GetCurrentCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void GetCurrentCommand::help(){
	try {
		m->mothurOut("The get.current command outputs the current files saved by mothur.\n");
		m->mothurOut("The get.current command has one parameter: clear.\n");
		m->mothurOut("The clear paramter is used to indicate which file types you would like to clear values for, multiple types can be separated by dashes.\n");
		m->mothurOut("The get.current command should be in the following format: \n");
		m->mothurOut("get.current() or get.current(clear=fasta-name-accnos)\n");
		
	}
	catch(exception& e) {
		m->errorOut(e, "GetCurrentCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************
GetCurrentCommand::~GetCurrentCommand(){}
//**********************************************************************************************************************

int GetCurrentCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//user wants to clear a type
		if (types.size() != 0) {
			for (int i = 0; i < types.size(); i++) {
				
				if (m->control_pressed) { break; }
				
				//look for file types
				if (types[i] == "fasta") {
					m->setFastaFile("");
				}else if (types[i] == "qfile") {
					m->setQualFile("");
				}else if (types[i] == "phylip") {
					m->setPhylipFile("");
				}else if (types[i] == "column") {
					m->setColumnFile("");
				}else if (types[i] == "list") {
					m->setListFile("");
				}else if (types[i] == "rabund") {
					m->setRabundFile("");
				}else if (types[i] == "sabund") {
					m->setSabundFile("");
				}else if (types[i] == "name") {
					m->setNameFile("");
				}else if (types[i] == "group") {
					m->setGroupFile("");
				}else if (types[i] == "order") {
					m->setOrderFile("");
				}else if (types[i] == "ordergroup") {
					m->setOrderGroupFile("");
				}else if (types[i] == "tree") {
					m->setTreeFile("");
				}else if (types[i] == "shared") {
					m->setSharedFile("");
				}else if (types[i] == "relabund") {
					m->setRelAbundFile("");
				}else if (types[i] == "design") {
					m->setDesignFile("");
				}else if (types[i] == "sff") {
					m->setSFFFile("");
				}else if (types[i] == "oligos") {
					m->setOligosFile("");
				}else if (types[i] == "accnos") {
					m->setAccnosFile("");
				}else if (types[i] == "taxonomy") {
						m->setTaxonomyFile("");
				}else if (types[i] == "all") {
					m->clearCurrentFiles();
				}else {
					m->mothurOut("[ERROR]: mothur does not save a current file for " + types[i]); m->mothurOutEndLine();
				}
			}
		}
		
		m->mothurOutEndLine(); m->mothurOut("Current files saved by mothur:"); m->mothurOutEndLine();
		m->printCurrentFiles();
		
		return 0;	
	}
	
	catch(exception& e) {
		m->errorOut(e, "GetCurrentCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************



