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
			}

			//you are doing a list and group shared
			if ((listfile != "") && (groupfile != "")) { globaldata->setFormat("shared"); }
			
			//you have not given a file
			if ((listfile == "") && (sharedfile == "") && (rabundfile == "") && (sabundfile == "")) {
				cout << "You must enter either a listfile, rabundfile, sabundfile or a sharedfile with the read.otu command. " << endl; abort = true; 
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
			if ((line != "") && (label != "")) { cout << "You cannot use both the line and label parameters at the same time. " << endl; abort = true; }
			
			orderfile = validParameter.validFile(parameters, "order", true);
			if (orderfile == "not open") { abort = true; }	
			else if (orderfile == "not found") { orderfile = ""; }
			else {  globaldata->setOrderFile(orderfile);	}
			
			
			if (abort == false) {
				//gets whichever one of the above is set
				filename = globaldata->inputFileName;
				read = new ReadOTUFile(filename);
			}

		}

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadOtuCommand class Function ReadOtuCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadOtuCommand class function ReadOtuCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
//**********************************************************************************************************************

void ReadOtuCommand::help(){
	try {
		cout << "The read.otu command must be run before you execute a collect.single, rarefaction.single, summary.single, " << "\n";
		cout << "collect.shared, rarefaction.shared or summary.shared command.   Mothur will generate a .list, .rabund and .sabund upon completion of the cluster command " << "\n";
		cout << "or you may use your own. The read.otu command parameter options are list, rabund, sabund, shared, group, order, line and label." << "\n";
		cout << "The read.otu command can be used in two ways.  The first is to read a list, rabund or sabund and run the collect.single, rarefaction.single or summary.single." << "\n";
		cout << "For this use the read.otu command should be in the following format: read.otu(list=yourListFile, order=yourOrderFile, label=yourLabels)." << "\n";
		cout << "The list, rabund or sabund parameter is required, but you may only use one of them." << "\n";
		cout << "The line and label parameters are optional but you may not use both the line and label parameters at the same time." << "\n";
		cout << "The label and line parameters are used to read specific lines in your input." << "\n";
		cout << "The second way to use the read.otu command is to read a list and a group, or a shared so you can use the collect.shared, rarefaction.shared or summary.shared commands." << "\n";
		cout << "In this case the read.otu command should be in the following format: read.otu(list=yourListFile, group=yourGroupFile, line=yourLines) or read.otu(shared=yourSharedFile).  " << "\n";
		cout << "The list parameter and group paramaters or the shared paremeter is required. When using the command the second way with a list and group file read.otu command parses the .list file" << "\n";
		cout << "and separates it into groups.  It outputs a .shared file containing the OTU information for each group. The read.otu command also outputs a .list file for each group. " << "\n";
		cout << "Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListfile)." << "\n" << "\n";

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadOtuCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadOtuCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}



//**********************************************************************************************************************

ReadOtuCommand::~ReadOtuCommand(){
	if (abort == false) {  delete read;  }
}

//**********************************************************************************************************************

int ReadOtuCommand::execute(){
	try {
	
		if (abort == true) {	return 0;	}
		
		read->read(&*globaldata); 
		if (globaldata->getFormat() == "shared") {
			groupMap->readMap();
			
			//if (globaldata->gGroupmap != NULL) { delete globaldata->gGroupmap;  }
			globaldata->gGroupmap = groupMap;		
			shared = new SharedCommand();
			shared->execute();

			parselist = new ParseListCommand();
			parselist->execute();
			
			//change format to shared  to speed up commands
			globaldata->setFormat("sharedfile");
			globaldata->setListFile("");
			globaldata->setGroupFile("");
			globaldata->setSharedFile(getRootName(filename) + "shared");
		}
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the ReadOtuCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the ReadOtuCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}
//**********************************************************************************************************************
