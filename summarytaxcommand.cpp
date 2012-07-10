/*
 *  summarytaxcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 9/23/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "summarytaxcommand.h"
#include "phylosummary.h"

//**********************************************************************************************************************
vector<string> SummaryTaxCommand::setParameters(){	
	try {
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(ptaxonomy);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pname);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pgroup);
		CommandParameter preftaxonomy("reftaxonomy", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(preftaxonomy);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryTaxCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SummaryTaxCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The summary.tax command reads a taxonomy file and an optional name file, and summarizes the taxonomy information.\n";
		helpString += "The summary.tax command parameters are taxonomy, group and name. taxonomy is required, unless you have a valid current taxonomy file.\n";
		helpString += "The name parameter allows you to enter a name file associated with your taxonomy file. \n";
		helpString += "The group parameter allows you add a group file so you can have the summary totals broken up by group.\n";
		helpString += "The reftaxonomy parameter allows you give the name of the reference taxonomy file used when you classified your sequences. It is not required, but providing it will keep the rankIDs in the summary file static.\n";
		helpString += "The summary.tax command should be in the following format: \n";
		helpString += "summary.tax(taxonomy=yourTaxonomyFile) \n";
		helpString += "Note: No spaces between parameter labels (i.e. taxonomy), '=' and parameters (i.e.yourTaxonomyFile).\n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryTaxCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string SummaryTaxCommand::getOutputFileNameTag(string type, string inputName=""){	
	try {
        string outputFileName = "";
		map<string, vector<string> >::iterator it;
        
        //is this a type this command creates
        it = outputTypes.find(type);
        if (it == outputTypes.end()) {  m->mothurOut("[ERROR]: this command doesn't create a " + type + " output file.\n"); }
        else {
            if (type == "summary")            {   outputFileName =  "tax.summary";   }
            else { m->mothurOut("[ERROR]: No definition for type " + type + " output file tag.\n"); m->control_pressed = true;  }
        }
        return outputFileName;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryTaxCommand", "getOutputFileNameTag");
		exit(1);
	}
}
//**********************************************************************************************************************
SummaryTaxCommand::SummaryTaxCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["summary"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryTaxCommand", "SummaryTaxCommand");
		exit(1);
	}
}
//***************************************************************************************************************

SummaryTaxCommand::SummaryTaxCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
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
				it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
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
				
				it = parameters.find("reftaxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["reftaxonomy"] = inputDir + it->second;		}
				}
				
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["summary"] = tempOutNames;
			
			//check for required parameters
			taxfile = validParameter.validFile(parameters, "taxonomy", true);
			if (taxfile == "not open") { abort = true; }
			else if (taxfile == "not found") { 				
				taxfile = m->getTaxonomyFile(); 
				if (taxfile != "") { m->mothurOut("Using " + taxfile + " as input file for the taxonomy parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current taxonomy file and the taxonomy parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setTaxonomyFile(taxfile); }	
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") { namefile = "";  }	
			else { m->setNameFile(namefile); }
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { groupfile = ""; abort = true; }
			else if (groupfile == "not found") { groupfile = ""; }
			else { m->setGroupFile(groupfile); }
			
			refTaxonomy = validParameter.validFile(parameters, "reftaxonomy", true);
			if (refTaxonomy == "not found") { refTaxonomy = ""; m->mothurOut("reftaxonomy is not required, but if given will keep the rankIDs in the summary file static."); m->mothurOutEndLine(); }
			else if (refTaxonomy == "not open") { refTaxonomy = ""; abort = true; }
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(taxfile); //if user entered a file with a path then preserve it	
			}
			
			if (namefile == "") {
				vector<string> files; files.push_back(taxfile);
				parser.getNameFile(files);
			}
			
		}
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryTaxCommand", "SummaryTaxCommand");
		exit(1);
	}
}
//***************************************************************************************************************

int SummaryTaxCommand::execute(){
	try{
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		int start = time(NULL);
		
		PhyloSummary* taxaSum;
		if (refTaxonomy != "") {
			taxaSum = new PhyloSummary(refTaxonomy, groupfile);
		}else {
			taxaSum = new PhyloSummary(groupfile);
		}
		
		if (m->control_pressed) { delete taxaSum; return 0; }
		
		int numSeqs = 0;
		if (namefile == "") { numSeqs = taxaSum->summarize(taxfile);  }
		else {
			map<string, vector<string> > nameMap;
			map<string, vector<string> >::iterator itNames;
			m->readNames(namefile, nameMap);
			
			if (m->control_pressed) { delete taxaSum; return 0; }
			
			ifstream in;
			m->openInputFile(taxfile, in);
			
			//read in users taxonomy file and add sequences to tree
			string name, taxon;
			
			while(!in.eof()){
				in >> name >> taxon; m->gobble(in);
				
				itNames = nameMap.find(name);
				
				if (itNames == nameMap.end()) { 
					m->mothurOut("[ERROR]: " + name + " is not in your name file please correct."); m->mothurOutEndLine(); exit(1);
				}else{
					for (int i = 0; i < itNames->second.size(); i++) { 
						numSeqs++;
						taxaSum->addSeqToTree(itNames->second[i], taxon);  //add it as many times as there are identical seqs
					}
					itNames->second.clear();
					nameMap.erase(itNames->first);
				}
			}
			in.close();
		}
		
		if (m->control_pressed) {  delete taxaSum; return 0; }
		
		//print summary file
		ofstream outTaxTree;
		string summaryFile = outputDir + m->getRootName(m->getSimpleName(taxfile)) + getOutputFileNameTag("summary");
		m->openOutputFile(summaryFile, outTaxTree);
		taxaSum->print(outTaxTree);
		outTaxTree.close();
		
		delete taxaSum;
		
		if (m->control_pressed) {  m->mothurRemove(summaryFile); return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("It took " + toString(time(NULL) - start) + " secs to create the summary file for " + toString(numSeqs) + " sequences."); m->mothurOutEndLine(); m->mothurOutEndLine();
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		m->mothurOut(summaryFile); m->mothurOutEndLine();	outputNames.push_back(summaryFile); outputTypes["summary"].push_back(summaryFile);
		m->mothurOutEndLine();
					
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SummaryTaxCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************/


