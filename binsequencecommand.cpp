/*
 *  binsequencecommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/3/09.
 *  Copyright 2009 Schloss Lab UMASS Amhers. All rights reserved.
 *
 */

#include "binsequencecommand.h"

//**********************************************************************************************************************
vector<string> BinSeqCommand::getValidParameters(){	
	try {
		string AlignArray[] =  {"fasta","label","name", "group","outputdir","inputdir"};
		vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "BinSeqCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> BinSeqCommand::getRequiredParameters(){	
	try {
		string AlignArray[] =  {"fasta"};
		vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "BinSeqCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> BinSeqCommand::getRequiredFiles(){	
	try {
		string AlignArray[] =  {"list"};
		vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "BinSeqCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************
BinSeqCommand::BinSeqCommand(){	
	try {
		abort = true; calledHelp = true; 
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "BinSeqCommand", "BinSeqCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
BinSeqCommand::BinSeqCommand(string option) {
	try {
		globaldata = GlobalData::getInstance();
		abort = false; calledHelp = false;   
		allLines = 1;
		labels.clear();
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			//valid paramters for this command
			string AlignArray[] =  {"fasta","label","name", "group","outputdir","inputdir"};
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
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
			outputTypes["fasta"] = tempOutNames;
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(globaldata->getListFile()); //if user entered a file with a path then preserve it	
			}

			
			//make sure the user has already run the read.otu command
			if (globaldata->getListFile() == "") { 
				m->mothurOut("You must read a listfile before running the bin.seqs command."); 
				m->mothurOutEndLine(); 
				abort = true; 
			}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
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
			}

			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not found") { m->mothurOut("fasta is a required parameter for the bin.seqs command.");  m->mothurOutEndLine(); abort = true; }
			else if (fastafile == "not open") { abort = true; }	
		
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			//if the user has not specified any labels use the ones from read.otu
			if (label == "") {  
				allLines = globaldata->allLines; 
				labels = globaldata->labels; 
			}
			
			namesfile = validParameter.validFile(parameters, "name", true);
			if (namesfile == "not open") { abort = true; }	
			else if (namesfile == "not found") { namesfile = ""; }

			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }
			else if (groupfile == "not found") { groupfile = ""; }
			
			if (abort == false) { 
//				m->openInputFile(fastafile, in);
				fasta = new FastaMap();
				if (groupfile != "") {
					groupMap = new GroupMap(groupfile);
					
					int error = groupMap->readMap();
					if (error == 1) { delete groupMap; abort = true; }
				}
			}
	
		}
	}
	catch(exception& e) {
		m->errorOut(e, "BinSeqCommand", "BinSeqCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void BinSeqCommand::help(){
	try {
		m->mothurOut("The bin.seqs command can only be executed after a successful read.otu command of a listfile.\n");
		m->mothurOut("The bin.seqs command parameters are fasta, name, label and group.  The fasta parameter is required.\n");
		m->mothurOut("The label parameter allows you to select what distance levels you would like a output files created for, and are separated by dashes.\n");
		m->mothurOut("The bin.seqs command should be in the following format: bin.seqs(fasta=yourFastaFile, name=yourNamesFile, group=yourGroupFile, label=yourLabels).\n");
		m->mothurOut("Example bin.seqs(fasta=amazon.fasta, group=amazon.groups, name=amazon.names).\n");
		m->mothurOut("The default value for label is all lines in your inputfile.\n");
		m->mothurOut("The bin.seqs command outputs a .fasta file for each distance you specify appending the OTU number to each name.\n");
		m->mothurOut("If you provide a groupfile, then it also appends the sequences group to the name.\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "BinSeqCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

BinSeqCommand::~BinSeqCommand(){
	//made new in execute
	if (abort == false) {
		delete input;  globaldata->ginput = NULL;
		delete read;
		globaldata->gListVector = NULL;
		delete fasta;
		if (groupfile != "") {  delete groupMap;  globaldata->gGroupmap = NULL; }
	}
}

//**********************************************************************************************************************

int BinSeqCommand::execute(){
	try {
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
	
		int error = 0;
		
		//read fastafile
		fasta->readFastaFile(fastafile);
		
		
		//set format to list so input can get listvector
//		globaldata->setFormat("list");
		
		//if user gave a namesfile then use it
		if (namesfile != "") {
			readNamesFile();
		}
		
		//read list file
		read = new ReadOTUFile(globaldata->getListFile());	
		read->read(&*globaldata); 
		
		input = globaldata->ginput;
		list = globaldata->gListVector;
		string lastLabel = list->getLabel();
		
		if (m->control_pressed) {  return 0; }
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;

				
		while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if(m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());		} return 0; }	
			
			if(allLines == 1 || labels.count(list->getLabel()) == 1){
				
				error = process(list);	
				if (error == 1) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());		} return 0; }	
							
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
			}
			
			if ((m->anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = list->getLabel();
				
				delete list;
				list = input->getListVector(lastLabel);
				
				error = process(list);	
				if (error == 1) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());		} return 0; }
													
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
				
				//restore real lastlabel to save below
				list->setLabel(saveLabel);
			}
			
			lastLabel = list->getLabel();			
			
			delete list;
			list = input->getListVector();
		}
		
		if(m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());		} return 0; }	

		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			m->mothurOut("Your file does not include the label " + *it); 
			if (processedLabels.count(lastLabel) != 1) {
				m->mothurOut(". I will use " + lastLabel + "."); m->mothurOutEndLine();
				needToRun = true;
			}else {
				m->mothurOut(". Please refer to " + lastLabel + ".");  m->mothurOutEndLine();
			}
		}
		
		//run last label if you need to
		if (needToRun == true)  {
			if (list != NULL) {		delete list;	}
			list = input->getListVector(lastLabel);
				
			error = process(list);	
			if (error == 1) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());		} return 0; }
			
			delete list;  
		}
		
		if(m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());		} return 0; }	

		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();

		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "BinSeqCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
void BinSeqCommand::readNamesFile() {
	try {
		vector<string> dupNames;
		m->openInputFile(namesfile, inNames);
		
		string name, names, sequence;
	
		while(inNames){
			inNames >> name;			//read from first column  A
			inNames >> names;		//read from second column  A,B,C,D
			
			dupNames.clear();
			
			//parse names into vector
			m->splitAtComma(names, dupNames);
			
			//store names in fasta map
			sequence = fasta->getSequence(name);
			for (int i = 0; i < dupNames.size(); i++) {
				fasta->push_back(dupNames[i], sequence);
			}
		
			m->gobble(inNames);
		}
		inNames.close();

	}
	catch(exception& e) {
		m->errorOut(e, "BinSeqCommand", "readNamesFile");
		exit(1);
	}
}
//**********************************************************************************************************************
//return 1 if error, 0 otherwise
int BinSeqCommand::process(ListVector* list) {
	try {
				string binnames, name, sequence;
				
				string outputFileName = outputDir + m->getRootName(m->getSimpleName(globaldata->getListFile())) + list->getLabel() + ".fasta";
				m->openOutputFile(outputFileName, out);
				
				//save to output list of output file names
				outputNames.push_back(outputFileName);  outputTypes["fasta"].push_back(outputFileName);

				m->mothurOut(list->getLabel()); m->mothurOutEndLine();
				
				//for each bin in the list vector
				for (int i = 0; i < list->size(); i++) {
					
					if (m->control_pressed) {  return 1; }
					
					binnames = list->get(i);
					while (binnames.find_first_of(',') != -1) { 
						name = binnames.substr(0,binnames.find_first_of(','));
						binnames = binnames.substr(binnames.find_first_of(',')+1, binnames.length());
						
						//do work for that name
						sequence = fasta->getSequence(name);
						if (sequence != "not found") {
							//if you don't have groups
							if (groupfile == "") {
								name = name + "\t" + toString(i+1);
								out << ">" << name << endl;
								out << sequence << endl;
							}else {//if you do have groups
								string group = groupMap->getGroup(name);
								if (group == "not found") {  
									m->mothurOut(name + " is missing from your group file. Please correct. ");  m->mothurOutEndLine();
									return 1;
								}else{
									name = name + "\t" + group + "\t" + toString(i+1);
									out << ">" << name << endl;
									out << sequence << endl;
								}
							}
						}else { 
							m->mothurOut(name + " is missing from your fasta or name file. Please correct. "); m->mothurOutEndLine();
							return 1;
						}
						
					}
					
					//get last name
					sequence = fasta->getSequence(binnames);
					if (sequence != "not found") {
						//if you don't have groups
						if (groupfile == "") {
							binnames = binnames + "\t" + toString(i+1);
							out << ">" << binnames << endl;
							out << sequence << endl;
						}else {//if you do have groups
							string group = groupMap->getGroup(binnames);
							if (group == "not found") {  
								m->mothurOut(binnames + " is missing from your group file. Please correct. "); m->mothurOutEndLine();
								return 1;
							}else{
								binnames = binnames + "\t" + group + "\t" + toString(i+1);
								out << ">" << binnames << endl;
								out << sequence << endl;
							}
						}
					}else { 
						m->mothurOut(binnames + " is missing from your fasta or name file. Please correct. "); m->mothurOutEndLine();
						return 1;
					}
				}
					
				out.close();
				return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "BinSeqCommand", "process");
		exit(1);
	}
}
//**********************************************************************************************************************


