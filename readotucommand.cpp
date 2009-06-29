/*
 *  readotu.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 1/20/09.
 *  Copyright 2009 Schloss Lab UMASS AMherst. All rights reserved.
 *
 */

#include "readotucommand.h"

//**********************************************************************************************************************
ReadOtuCommand::ReadOtuCommand(string option){
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		allLines = 1;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"list","order","shared", "line", "label","group","sabund", "rabund"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string, string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			globaldata->newRead();
			
			//check for required parameters
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { abort = true; }
			else if (listfile == "not found") { listfile = ""; }	
			else {  globaldata->setListFile(listfile);  globaldata->setFormat("list"); 	}
			
			sabundfile = validParameter.validFile(parameters, "sabund", true);
			if (sabundfile == "not open") { abort = true; }	
			else if (sabundfile == "not found") { sabundfile = ""; }
			else {  globaldata->setSabundFile(sabundfile); globaldata->setFormat("sabund");	}

			rabundfile = validParameter.validFile(parameters, "rabund", true);
			if (rabundfile == "not open") { abort = true; }	
			else if (rabundfile == "not found") { rabundfile = ""; }
			else {  globaldata->setRabundFile(rabundfile);	globaldata->setFormat("rabund");}
			
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { abort = true; }	
			else if (sharedfile == "not found") { sharedfile = ""; }
			else {  globaldata->setSharedFile(sharedfile); globaldata->setFormat("sharedfile");	}
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }	
			else if (groupfile == "not found") { groupfile = ""; }
			else {  
				globaldata->setGroupFile(groupfile); 
				groupMap = new GroupMap(groupfile);
				groupMap->readMap();
				globaldata->gGroupmap = groupMap;
			}

			//you are doing a list and group shared
			if ((listfile != "") && (groupfile != "")) { globaldata->setFormat("shared"); }
			
			//you have not given a file
			if ((listfile == "") && (sharedfile == "") && (rabundfile == "") && (sabundfile == "")) {
				mothurOut("You must enter either a listfile, rabundfile, sabundfile or a sharedfile with the read.otu command. "); mothurOutEndLine(); abort = true; 
			}
		
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			line = validParameter.validFile(parameters, "line", false);				
			if (line == "not found") { line = ""; }
			else { 
				if(line != "all") {  splitAtDash(line, lines);  allLines = 0;  }
				else { allLines = 1;  }
				globaldata->lines = lines;
			}
			
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
				globaldata->labels = labels;
			}
			
			globaldata->allLines = allLines;
			
			//make sure user did not use both the line and label parameters
			if ((line != "") && (label != "")) { mothurOut("You cannot use both the line and label parameters at the same time. "); mothurOutEndLine(); abort = true; }
			
			orderfile = validParameter.validFile(parameters, "order", true);
			if (orderfile == "not open") { abort = true; }	
			else if (orderfile == "not found") { orderfile = ""; }
			else {  globaldata->setOrderFile(orderfile);	}
			
			
			if (abort == false) {
				//gets whichever one of the above is set
				filename = globaldata->inputFileName;
			}

		}

	}
	catch(exception& e) {
		errorOut(e, "ReadOtuCommand", "ReadOtuCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void ReadOtuCommand::help(){
	try {
		mothurOut("The read.otu command must be run before you execute a collect.single, rarefaction.single, summary.single, \n");
		mothurOut("collect.shared, rarefaction.shared or summary.shared command.   Mothur will generate a .list, .rabund and .sabund upon completion of the cluster command \n");
		mothurOut("or you may use your own. The read.otu command parameter options are list, rabund, sabund, shared, group, order, line and label.\n");
		mothurOut("The read.otu command can be used in two ways.  The first is to read a list, rabund or sabund and run the collect.single, rarefaction.single or summary.single.\n");
		mothurOut("For this use the read.otu command should be in the following format: read.otu(list=yourListFile, order=yourOrderFile, label=yourLabels).\n");
		mothurOut("The list, rabund or sabund parameter is required, but you may only use one of them.\n");
		mothurOut("The line and label parameters are optional but you may not use both the line and label parameters at the same time.\n");
		mothurOut("The label and line parameters are used to read specific lines in your input.\n");
		mothurOut("The second way to use the read.otu command is to read a list and a group, or a shared so you can use the collect.shared, rarefaction.shared or summary.shared commands.\n");
		mothurOut("In this case the read.otu command should be in the following format: read.otu(list=yourListFile, group=yourGroupFile, line=yourLines) or read.otu(shared=yourSharedFile).  \n");
		mothurOut("The list parameter and group paramaters or the shared paremeter is required. When using the command the second way with a list and group file read.otu command parses the .list file\n");
		mothurOut("and separates it into groups.  It outputs a .shared file containing the OTU information for each group. The read.otu command also outputs a .list file for each group. \n");
		mothurOut("Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListfile).\n\n");

	}
	catch(exception& e) {
		errorOut(e, "ReadOtuCommand", "help");
		exit(1);
	}
}



//**********************************************************************************************************************

ReadOtuCommand::~ReadOtuCommand(){
	}

//**********************************************************************************************************************

int ReadOtuCommand::execute(){
	try {
	
		if (abort == true) {	return 0;	}
		
		if (globaldata->getFormat() == "shared") {
			
			shared = new SharedCommand();
			shared->execute();
			delete shared;
				
			//change format to shared  to speed up commands
			globaldata->setFormat("sharedfile");
			globaldata->setListFile("");
			globaldata->setGroupFile("");
			globaldata->setSharedFile(getRootName(filename) + "shared");
		}
		return 0;
	}
	catch(exception& e) {
		errorOut(e, "ReadOtuCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
