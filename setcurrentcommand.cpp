/*
 *  setcurrentcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 3/16/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "setcurrentcommand.h"

//**********************************************************************************************************************
vector<string> SetCurrentCommand::getValidParameters(){	
	try {
		string Array[] =  {"outputdir","inputdir","clear", "phylip", "column", "list","rabund","sabund","name","group",
			"design","order","tree","shared","ordergroup","relabund","fasta","qfile","sff","oligos","accnos","taxonomy"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SetCurrentCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
SetCurrentCommand::SetCurrentCommand(){	
	try {
		abort = true; calledHelp = true; 
	}
	catch(exception& e) {
		m->errorOut(e, "SetCurrentCommand", "SetCurrentCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> SetCurrentCommand::getRequiredParameters(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SetCurrentCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> SetCurrentCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SetCurrentCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************
SetCurrentCommand::SetCurrentCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"outputdir","inputdir","clear", "phylip", "column", "list","rabund","sabund","name","group",
				"design","order","tree","shared","ordergroup","relabund","fasta","qfile","sff","oligos","accnos","taxonomy"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("phylip");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["phylip"] = inputDir + it->second;		}
				}
				
				it = parameters.find("column");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["column"] = inputDir + it->second;		}
				}
				
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
				
				it = parameters.find("rabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["rabund"] = inputDir + it->second;		}
				}
				
				it = parameters.find("sabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["sabund"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
				
				it = parameters.find("design");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["design"] = inputDir + it->second;		}
				}
				
				it = parameters.find("order");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["order"] = inputDir + it->second;		}
				}
				
				it = parameters.find("tree");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["tree"] = inputDir + it->second;		}
				}
				
				it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
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
				
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("qfile");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["qfile"] = inputDir + it->second;		}
				}
				
				it = parameters.find("sff");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["sff"] = inputDir + it->second;		}
				}
				
				it = parameters.find("oligos");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["oligos"] = inputDir + it->second;		}
				}
				
				it = parameters.find("accnos");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["accnos"] = inputDir + it->second;		}
				}

				it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
				}
			}
			
			//check for parameters
			phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not open") { m->mothurOut("Ignoring: " + parameters["phylip"]); m->mothurOutEndLine(); phylipfile = ""; }
			else if (phylipfile == "not found") {  phylipfile = "";  }	
			if (phylipfile != "") { m->setPhylipFile(phylipfile); }
			
			columnfile = validParameter.validFile(parameters, "column", true);
			if (columnfile == "not open") { m->mothurOut("Ignoring: " + parameters["column"]); m->mothurOutEndLine(); columnfile = ""; }
			else if (columnfile == "not found") {  columnfile = "";  }	
			if (columnfile != "") { m->setColumnFile(columnfile); }
			
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { m->mothurOut("Ignoring: " + parameters["list"]); m->mothurOutEndLine(); listfile = ""; }
			else if (listfile == "not found") {  listfile = "";  }	
			if (listfile != "") { m->setListFile(listfile); }
			
			rabundfile = validParameter.validFile(parameters, "rabund", true);
			if (rabundfile == "not open") { m->mothurOut("Ignoring: " + parameters["rabund"]); m->mothurOutEndLine(); rabundfile = ""; }
			else if (rabundfile == "not found") {  rabundfile = "";  }	
			if (rabundfile != "") { m->setRabundFile(rabundfile); }
			
			sabundfile = validParameter.validFile(parameters, "sabund", true);
			if (sabundfile == "not open") { m->mothurOut("Ignoring: " + parameters["sabund"]); m->mothurOutEndLine(); sabundfile = ""; }
			else if (sabundfile == "not found") {  sabundfile = "";  }	
			if (sabundfile != "") { m->setSabundFile(sabundfile); }
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { m->mothurOut("Ignoring: " + parameters["name"]); m->mothurOutEndLine(); namefile = ""; }
			else if (namefile == "not found") {  namefile = "";  }	
			if (namefile != "") { m->setNameFile(namefile); }
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { m->mothurOut("Ignoring: " + parameters["group"]); m->mothurOutEndLine(); groupfile = ""; }
			else if (groupfile == "not found") {  groupfile = "";  }
			if (groupfile != "") { m->setGroupFile(groupfile); }
			
			designfile = validParameter.validFile(parameters, "design", true);
			if (designfile == "not open") { m->mothurOut("Ignoring: " + parameters["design"]); m->mothurOutEndLine(); designfile = ""; }
			else if (designfile == "not found") {  designfile = "";  }	
			if (designfile != "") { m->setDesignFile(designfile); }
			
			orderfile = validParameter.validFile(parameters, "order", true);
			if (orderfile == "not open") { m->mothurOut("Ignoring: " + parameters["order"]); m->mothurOutEndLine(); orderfile = ""; }
			else if (orderfile == "not found") {  orderfile = "";  }
			if (orderfile != "") { m->setOrderFile(orderfile); }
			
			treefile = validParameter.validFile(parameters, "tree", true);
			if (treefile == "not open") { m->mothurOut("Ignoring: " + parameters["tree"]); m->mothurOutEndLine(); treefile = ""; }
			else if (treefile == "not found") {  treefile = "";  }	
			if (treefile != "") { m->setTreeFile(treefile); }
			
			sharedfile = validParameter.validFile(parameters, "shared", true);
			if (sharedfile == "not open") { m->mothurOut("Ignoring: " + parameters["shared"]); m->mothurOutEndLine(); sharedfile = ""; }
			else if (sharedfile == "not found") {  sharedfile = "";  }	
			if (sharedfile != "") { m->setSharedFile(sharedfile); }
			
			ordergroupfile = validParameter.validFile(parameters, "ordergroup", true);
			if (ordergroupfile == "not open") { m->mothurOut("Ignoring: " + parameters["ordergroup"]); m->mothurOutEndLine(); ordergroupfile = ""; }
			else if (ordergroupfile == "not found") {  ordergroupfile = "";  }	
			if (ordergroupfile != "") { m->setOrderGroupFile(ordergroupfile); }
			
			relabundfile = validParameter.validFile(parameters, "relabund", true);
			if (relabundfile == "not open") { m->mothurOut("Ignoring: " + parameters["relabund"]); m->mothurOutEndLine(); relabundfile = ""; }
			else if (relabundfile == "not found") {  relabundfile = "";  }	
			if (relabundfile != "") { m->setRelAbundFile(relabundfile); }
			
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { m->mothurOut("Ignoring: " + parameters["fasta"]); m->mothurOutEndLine(); fastafile = ""; }
			else if (fastafile == "not found") {  fastafile = "";  }
			if (fastafile != "") { m->setFastaFile(fastafile); }
			
			qualfile = validParameter.validFile(parameters, "qfile", true);
			if (qualfile == "not open") { m->mothurOut("Ignoring: " + parameters["qfile"]); m->mothurOutEndLine(); qualfile = ""; }
			else if (qualfile == "not found") {  qualfile = "";  }	
			if (qualfile != "") { m->setQualFile(qualfile); }
			
			sfffile = validParameter.validFile(parameters, "sff", true);
			if (sfffile == "not open") { m->mothurOut("Ignoring: " + parameters["sff"]); m->mothurOutEndLine(); sfffile = ""; }
			else if (sfffile == "not found") {  sfffile = "";  }	
			if (sfffile != "") { m->setSFFFile(sfffile); }
			
			oligosfile = validParameter.validFile(parameters, "oligos", true);
			if (oligosfile == "not open") { m->mothurOut("Ignoring: " + parameters["oligos"]); m->mothurOutEndLine(); oligosfile = ""; }
			else if (oligosfile == "not found") {  oligosfile = "";  }	
			if (oligosfile != "") { m->setOligosFile(oligosfile); }
			
			accnosfile = validParameter.validFile(parameters, "accnos", true);
			if (accnosfile == "not open") { m->mothurOut("Ignoring: " + parameters["accnos"]); m->mothurOutEndLine(); accnosfile = ""; }
			else if (accnosfile == "not found") {  accnosfile = "";  }	
			if (accnosfile != "") { m->setAccnosFile(accnosfile); }
			
			taxonomyfile = validParameter.validFile(parameters, "taxonomy", true);
			if (taxonomyfile == "not open") { m->mothurOut("Ignoring: " + parameters["taxonomy"]); m->mothurOutEndLine(); taxonomyfile = ""; }
			else if (taxonomyfile == "not found") {  taxonomyfile = "";  }	
			if (taxonomyfile != "") { m->setTaxonomyFile(taxonomyfile); }
			
			clearTypes = validParameter.validFile(parameters, "clear", false);			
			if (clearTypes == "not found") { clearTypes = ""; }
			else { m->splitAtDash(clearTypes, types);	}
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "SetCurrentCommand", "SetCurrentCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void SetCurrentCommand::help(){
	try {
		m->mothurOut("The set.current command allows you to set the current files saved by mothur.\n");
		m->mothurOut("The set.current command parameters are: clear, phylip, column, list, rabund, sabund, name, group, design, order, tree, shared, ordergroup, relabund, fasta, qfile, sff, oligos, accnos, taxonomy.\n");
		m->mothurOut("The clear paramter is used to indicate which file types you would like to clear values for, multiple types can be separated by dashes.\n");
		m->mothurOut("The set.current command should be in the following format: \n");
		m->mothurOut("set.current(fasta=yourFastaFile) or set.current(fasta=amazon.fasta, clear=name-accnos)\n");
		
	}
	catch(exception& e) {
		m->errorOut(e, "SetCurrentCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************
SetCurrentCommand::~SetCurrentCommand(){}
//**********************************************************************************************************************

int SetCurrentCommand::execute(){
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
		m->errorOut(e, "SetCurrentCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************



