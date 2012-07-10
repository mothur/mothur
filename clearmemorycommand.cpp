/*
 *  clearmemorycommand.cpp
 *  Mothur
 *
 *  Created by westcott on 7/6/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "clearmemorycommand.h"
#include "referencedb.h"

//**********************************************************************************************************************
vector<string> ClearMemoryCommand::setParameters(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ClearMemoryCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ClearMemoryCommand::getOutputFileNameTag(string type, string inputName=""){	
	try {
        string outputFileName = "";
		map<string, vector<string> >::iterator it;
        
        //is this a type this command creates
        it = outputTypes.find(type);
        if (it == outputTypes.end()) {  m->mothurOut("[ERROR]: this command doesn't create a " + type + " output file.\n"); }
        else {
            m->mothurOut("[ERROR]: No definition for type " + type + " output file tag.\n"); m->control_pressed = true; 
        }
        return outputFileName;
	}
	catch(exception& e) {
		m->errorOut(e, "ClearMemoryCommand", "getOutputFileNameTag");
		exit(1);
	}
}
//**********************************************************************************************************************
string ClearMemoryCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The clear.memory command removes saved reference data from memory.\n";
		helpString += "The clear.memory command should be in the following format: clear.memory().\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ClearMemoryCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************

ClearMemoryCommand::ClearMemoryCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
	}
	catch(exception& e) {
		m->errorOut(e, "ClearMemoryCommand", "ClearMemoryCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

int ClearMemoryCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		ReferenceDB* rdb = ReferenceDB::getInstance();
		rdb->clearMemory();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ClearMemoryCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************/
