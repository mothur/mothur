/*
 *  GetlineCommand.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 1/30/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "getlinecommand.h"

//**********************************************************************************************************************
GetlineCommand::GetlineCommand(string option){
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			if (option != "") { mothurOut("There are no valid parameters for the get.line command."); mothurOutEndLine(); abort = true; }
			
			if ((globaldata->getListFile() == "") && (globaldata->getRabundFile() == "") && (globaldata->getSabundFile() == "")) { mothurOut("You must read a list, sabund or rabund before you can use the get.line command."); mothurOutEndLine(); abort = true; }				
		}

	}
	catch(exception& e) {
		errorOut(e, "GetlineCommand", "GetlineCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void GetlineCommand::help(){
	try {
		mothurOut("The get.line command can only be executed after a successful read.otu command.\n");
		mothurOut("You may not use any parameters with the get.line command.\n");
		mothurOut("The get.line command should be in the following format: \n");
		mothurOut("get.line()\n");
		mothurOut("Example get.line().\n");
	}
	catch(exception& e) {
		errorOut(e, "GetlineCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

GetlineCommand::~GetlineCommand(){
}

//**********************************************************************************************************************

int GetlineCommand::execute(){
	try {
	
		if (abort == true) { return 0; }
	
		ifstream in;
		openInputFile(globaldata->inputFileName, in);
		string label;
		int numBins = 0;
		int count = -1;
		int line = 1;
		while(in.good()) {
			if(count > numBins)
				count = 0;
			if(count == 0) {
				mothurOut(toString(line)); mothurOutEndLine();
				in >> numBins;
				line++;
			}
			in >> label;
			count++;
		}
		
		in.close();
		return 0;		
	}

	catch(exception& e) {
		errorOut(e, "GetlineCommand", "execute");
		exit(1);
	}
}



