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
ReadOtuCommand::ReadOtuCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["rabund"] = tempOutNames;
		outputTypes["shared"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ReadOtuCommand", "ReadOtuCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
ReadOtuCommand::ReadOtuCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		allLines = 1;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			/*
			//valid paramters for this command
			string Array[] =  {"list","order","shared","relabund","label","group","sabund", "rabund","groups","ordergroup","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string, string>::iterator it;
		
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["rabund"] = tempOutNames;
			outputTypes["shared"] = tempOutNames;
			
			globaldata->newRead();
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
				
				it = parameters.find("order");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["order"] = inputDir + it->second;		}
				}
				
				it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
				
				it = parameters.find("sabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["sabund"] = inputDir + it->second;		}
				}
				
				it = parameters.find("rabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["rabund"] = inputDir + it->second;		}
				}
				
				it = parameters.find("ordergroup");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["ordergroup"] = inputDir + it->second;		}
				}
				
				it = parameters.find("relabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["relabund"] = inputDir + it->second;		}
				}
			}

			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}

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
			
			ordergroupfile = validParameter.validFile(parameters, "ordergroup", true);
			if (ordergroupfile == "not open") { abort = true; }	
			else if (ordergroupfile == "not found") { ordergroupfile = ""; }
			else {  globaldata->setOrderGroupFile(ordergroupfile);  }
			
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { abort = true; }	
			else if (sharedfile == "not found") { sharedfile = ""; }
			else {  globaldata->setSharedFile(sharedfile); globaldata->setFormat("sharedfile");	}
			
			relAbundfile = validParameter.validFile(parameters, "relabund", true);
			if (relAbundfile == "not open") { abort = true; }	
			else if (relAbundfile == "not found") { relAbundfile = ""; }
			else {  globaldata->setRelAbundFile(relAbundfile); globaldata->setFormat("relabund");	}

			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }	
			else if (groupfile == "not found") { groupfile = ""; }
			else {  
				globaldata->setGroupFile(groupfile); 
				groupMap = new GroupMap(groupfile);
				
				int error = groupMap->readMap();
				if (error == 1) { abort = true; }
				
				globaldata->gGroupmap = groupMap;
			}
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				m->splitAtDash(groups, Groups);
				globaldata->Groups = Groups;
			}

			//you are doing a list and group shared
			if ((listfile != "") && (groupfile != "")) { globaldata->setFormat("shared"); }
			
			//you have not given a file
			if ((listfile == "") && (sharedfile == "") && (rabundfile == "") && (sabundfile == "") && (relAbundfile == "")) {
				m->mothurOut("You must enter either a listfile, rabundfile, sabundfile, relabund or a sharedfile with the read.otu command. "); m->mothurOutEndLine(); abort = true; 
			}
		
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
				globaldata->labels = labels;
			}
			
			globaldata->allLines = allLines;
		
			orderfile = validParameter.validFile(parameters, "order", true);
			if (orderfile == "not open") { abort = true; }	
			else if (orderfile == "not found") { orderfile = ""; }
			else {  globaldata->setOrderFile(orderfile);	}
			
				
			if (abort == false) {
				//gets whichever one of the above is set
				filename = globaldata->inputFileName;
			}
			 */
		}

	}
	catch(exception& e) {
		m->errorOut(e, "ReadOtuCommand", "ReadOtuCommand");
		exit(1);
	}
}
///**********************************************************************************************************************

int ReadOtuCommand::execute(){
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		m->mothurOut(getHelpString()); m->mothurOutEndLine();
	
		/*
		if (globaldata->getFormat() == "shared") {
			
			shared = new SharedCommand(outputDir);
			int okay = shared->execute();
			
			//problem with shared
			if (okay == 1) {
				globaldata->setListFile("");
				globaldata->setGroupFile("");
				globaldata->setSharedFile("");
			}else { //shared command outputs the filenames
				//m->mothurOutEndLine();
				//m->mothurOut("Output File Name: "); m->mothurOutEndLine();
				//m->mothurOut(globaldata->getSharedFile()); m->mothurOutEndLine();	
				//m->mothurOutEndLine();
			}
			
			outputTypes = shared->getOutputFiles();
			
			//set rabund file as new current rabundfile
			string current = "";
			itTypes = outputTypes.find("rabund");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setRabundFile(current); }
			}
			
			itTypes = outputTypes.find("shared");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setSharedFile(current); }
			}			
			
			delete shared;
		}
		*/
				
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ReadOtuCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
