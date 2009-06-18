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
BinSeqCommand::BinSeqCommand(string option){
	try {
		globaldata = GlobalData::getInstance();
		abort = false;
		allLines = 1;
		lines.clear();
		labels.clear();
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string AlignArray[] =  {"fasta","line","label","name", "group"};
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string, string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//make sure the user has already run the read.otu command
			if (globaldata->getListFile() == "") { cout << "You must read a listfile before running the bin.seqs command." << endl; abort = true; }
			
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not found") { cout << "fasta is a required parameter for the bin.seqs command." << endl; abort = true; }
			else if (fastafile == "not open") { abort = true; }	
		
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			line = validParameter.validFile(parameters, "line", false);				
			if (line == "not found") { line = "";  }
			else { 
				if(line != "all") {  splitAtDash(line, lines);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			//make sure user did not use both the line and label parameters
			if ((line != "") && (label != "")) { cout << "You cannot use both the line and label parameters at the same time. " << endl; abort = true; }
			//if the user has not specified any line or labels use the ones from read.otu
			else if ((line == "") && (label == "")) {  
				allLines = globaldata->allLines; 
				labels = globaldata->labels; 
				lines = globaldata->lines;
			}
			
			namesfile = validParameter.validFile(parameters, "name", true);
			if (namesfile == "not open") { abort = true; }	
			else if (namesfile == "not found") { namesfile = ""; }

			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }
			else if (groupfile == "not found") { groupfile = ""; }
			
			if (abort == false) { 
				openInputFile(fastafile, in);
				fasta = new FastaMap();
				if (groupfile != "") {
					groupMap = new GroupMap(groupfile);
					groupMap->readMap();
				}
			}
	
		}
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the BinSeqCommand class Function BinSeqCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the BinSeqCommand class function BinSeqCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************

void BinSeqCommand::help(){
	try {
		cout << "The bin.seqs command can only be executed after a successful read.otu command of a listfile." << "\n";
		cout << "The bin.seqs command parameters are fasta, name, line, label and group.  The fasta parameter is required, and you may not use line and label at the same time." << "\n";
		cout << "The line and label allow you to select what distance levels you would like a output files created for, and are separated by dashes." << "\n";
		cout << "The bin.seqs command should be in the following format: bin.seqs(fasta=yourFastaFile, name=yourNamesFile, group=yourGroupFile, line=yourLines, label=yourLabels)." << "\n";
		cout << "Example bin.seqs(fasta=amazon.fasta, group=amazon.groups, line=1-3-5, name=amazon.names)." << "\n";
		cout << "The default value for line and label are all lines in your inputfile." << "\n";
		cout << "The bin.seqs command outputs a .fasta file for each distance you specify appending the OTU number to each name." << "\n";
		cout << "If you provide a groupfile, then it also appends the sequences group to the name." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile)." << "\n" << "\n";
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the BinSeqCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the BinSeqCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
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
		if (abort == true) {	return 0;	}
	
		int count = 1;
		int error = 0;
		
		//read fastafile
		fasta->readFastaFile(in);
		
		in.close();
		
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
		ListVector* lastList = list;
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		set<int> userLines = lines;

				
		while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {
			
			if(allLines == 1 || lines.count(count) == 1 || labels.count(list->getLabel()) == 1){
				
				error = process(list, count);	
				if (error == 1) { return 0; }	
							
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
				userLines.erase(count);
			}
			
			if ((anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastList->getLabel()) != 1)) {
				
				error = process(lastList, count);	
				if (error == 1) { return 0; }
													
				processedLabels.insert(lastList->getLabel());
				userLabels.erase(lastList->getLabel());
				
			}
			
			if (count != 1) { delete lastList; }
			lastList = list;			

			list = input->getListVector();
			count++;
		}
		
		
		//output error messages about any remaining user labels
		set<string>::iterator it;
		bool needToRun = false;
		for (it = userLabels.begin(); it != userLabels.end(); it++) {  
			cout << "Your file does not include the label "<< *it; 
			if (processedLabels.count(lastList->getLabel()) != 1) {
				cout << ". I will use " << lastList->getLabel() << "." << endl;
				needToRun = true;
			}else {
				cout << ". Please refer to " << lastList->getLabel() << "." << endl;
			}
		}
		
		//run last line if you need to
		if (needToRun == true)  {
			error = process(lastList, count);	
			if (error == 1) { return 0; }			
		}
		
		delete lastList;
		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the BinSeqCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the BinSeqCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************
void BinSeqCommand::readNamesFile() {
	try {
		vector<string> dupNames;
		openInputFile(namesfile, inNames);
		
		string name, names, sequence;
	
		while(inNames){
			inNames >> name;			//read from first column  A
			inNames >> names;		//read from second column  A,B,C,D
			
			dupNames.clear();
			
			//parse names into vector
			splitAtComma(names, dupNames);
			
			//store names in fasta map
			sequence = fasta->getSequence(name);
			for (int i = 0; i < dupNames.size(); i++) {
				fasta->push_back(dupNames[i], sequence);
			}
		
			gobble(inNames);
		}
		inNames.close();

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the BinSeqCommand class Function readNamesFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the BinSeqCommand class function readNamesFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************
//return 1 if error, 0 otherwise
int BinSeqCommand::process(ListVector* list, int count) {
	try {
				string binnames, name, sequence;
				string outputFileName = getRootName(globaldata->getListFile()) + list->getLabel() + ".fasta";
				openOutputFile(outputFileName, out);

				cout << list->getLabel() << '\t' << count << endl;
				
				//for each bin in the list vector
				for (int i = 0; i < list->size(); i++) {

					binnames = list->get(i);
					while (binnames.find_first_of(',') != -1) { 
						name = binnames.substr(0,binnames.find_first_of(','));
						binnames = binnames.substr(binnames.find_first_of(',')+1, binnames.length());
						
						//do work for that name
						sequence = fasta->getSequence(name);
						if (sequence != "not found") {
							//if you don't have groups
							if (groupfile == "") {
								name = name + "|" + toString(i+1);
								out << ">" << name << endl;
								out << sequence << endl;
							}else {//if you do have groups
								string group = groupMap->getGroup(name);
								if (group == "not found") {  
									cout << name << " is missing from your group file. Please correct. " << endl;
									remove(outputFileName.c_str());
									return 1;
								}else{
									name = name + "|" + group + "|" + toString(i+1);
									out << ">" << name << endl;
									out << sequence << endl;
								}
							}
						}else { 
							cout << name << " is missing from your fasta or name file. Please correct. " << endl; 
							remove(outputFileName.c_str());
							return 1;
						}
						
					}
					
					//get last name
					sequence = fasta->getSequence(binnames);
					if (sequence != "not found") {
						//if you don't have groups
						if (groupfile == "") {
							binnames = binnames + "|" + toString(i+1);
							out << ">" << binnames << endl;
							out << sequence << endl;
						}else {//if you do have groups
							string group = groupMap->getGroup(binnames);
							if (group == "not found") {  
								cout << binnames << " is missing from your group file. Please correct. " << endl;
								remove(outputFileName.c_str());
								return 1;
							}else{
								binnames = binnames + "|" + group + "|" + toString(i+1);
								out << ">" << binnames << endl;
								out << sequence << endl;
							}
						}
					}else { 
						cout << binnames << " is missing from your fasta or name file. Please correct. " << endl; 
						remove(outputFileName.c_str());
						return 1;
					}
				}
					
				out.close();
				return 0;

	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the BinSeqCommand class Function readNamesFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the BinSeqCommand class function readNamesFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************


