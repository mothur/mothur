/*
 *  systemcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "systemcommand.h"

//**********************************************************************************************************************
vector<string> SystemCommand::setParameters(){	
	try {
		     
        CurrentFile* current; current = CurrentFile::getInstance();
        outputdir = current->getOutputDir();
        
        abort = false; calledHelp = false;
				
		vector<string> myArray;
		
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SystemCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
SystemCommand::SystemCommand(string option) : Command()  {
	try {

		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
            setParameters();
            string optionCopy = option;
            string parameter = "";
            
            //if command is used parameter=command optionCopy=cp ....
            if (optionCopy.find("command=") != string::npos) { util.splitAtEquals(parameter, optionCopy); }
			
			//ValidParameters validParameter;
			if (parameter != "command") { command = option; }
			else { command = optionCopy; } //command= removed
			
            if ((command == "")) { m->mothurOut("[ERROR]: You must enter a command to run.\n"); abort = true; }
		}	

	}
	catch(exception& e) {
		m->errorOut(e, "SystemCommand", "SystemCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

string SystemCommand::getHelpString(){
	try {
		string helpString = "";
		helpString += "The system command allows you to execute a system command from within mothur.\n";
		helpString += "The system has no parameters.\n";
		helpString += "The system command should be in the following format: system(yourCommand).\n";
		helpString += "Example system(clear).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SystemCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

int SystemCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
        if (outputdir == "") { outputdir = "./"; }
        string redirectFileName = outputdir + "commandScreen.output";
        
		//if command contains a redirect don't add the redirect
		bool usedRedirect = false;
		if ((command.find('>')) == string::npos) {
			command += " > " + redirectFileName + " 2>&1";
			usedRedirect = true;
		}
       
        //m->mothurOut("[DEBUG]: command = '" + command + "'\n");
        
        if (m->getDebug()) { m->mothurOut("[DEBUG]: command = '" + command + "'\n"); }
		
        system(command.c_str());
  
		if (usedRedirect) {
			ifstream in; util.openInputFile(redirectFileName, in, "no error");
			
			string output = "";
			while(char c = in.get()){
				if(in.eof())		{	break;			}
				else				{	output += c;	}
			}
			in.close();
            
			m->mothurOut(output); m->mothurOutEndLine();
			util.mothurRemove(redirectFileName);
		}
		
		return 0;		
	}

	catch(exception& e) {
		m->errorOut(e, "SystemCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
