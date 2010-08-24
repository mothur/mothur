/*
 *  GetlabelCommand.cpp
 *  Mothur
 *
 *  Created by Thomas Ryabin on 1/30/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "getlabelcommand.h"

//**********************************************************************************************************************

GetlabelCommand::GetlabelCommand(string option)  {
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			if ((globaldata->getListFile() == "") && (globaldata->getRabundFile() == "") && (globaldata->getSabundFile() == "")) { m->mothurOut("You must read a list, sabund or rabund before you can use the get.label command."); m->mothurOutEndLine(); abort = true; }				
		}

	}
	catch(exception& e) {
		m->errorOut(e, "GetlabelCommand", "GetlabelCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void GetlabelCommand::help(){
	try {
		m->mothurOut("The get.label command can only be executed after a successful read.otu command.\n");
		m->mothurOut("You may not use any parameters with the get.label command.\n");
		m->mothurOut("The get.label command should be in the following format: \n");
		m->mothurOut("get.label()\n");
		m->mothurOut("Example get.label().\n");
	}
	catch(exception& e) {
		m->errorOut(e, "GetlabelCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

GetlabelCommand::~GetlabelCommand(){
}

//**********************************************************************************************************************

int GetlabelCommand::execute(){
	try {
		
		if (abort == true) { return 0; }
		
		filename = globaldata->inputFileName;
		ifstream in;
		m->openInputFile(filename, in);
		string label;
		int numBins = 0;
		int count = -1;
		while(in.good()) {
			
			if (m->control_pressed) { in.close();  return 0; }
			
			if(count > numBins)
				count = 0;
			if(count == 0) {
				m->mothurOut(label); m->mothurOutEndLine();
				in >> numBins;
			}
			in >> label;
			count++;
		}	
		
		in.close();
		return 0;	
	}

	catch(exception& e) {
		m->errorOut(e, "GetlabelCommand", "execute");
		exit(1);
	}
}

