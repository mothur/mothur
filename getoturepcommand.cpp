/*
 *  getoturepcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/6/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "getoturepcommand.h"

//**********************************************************************************************************************
GetOTURepCommand::GetOTURepCommand(string option){
	try{
		globaldata = GlobalData::getInstance();
		abort = false;
		allLines = 1;
		lines.clear();
		labels.clear();
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta","list","line","label","name", "group"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
		
			//check to make sure all parameters are valid for command
			for (map<string, string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//make sure the user has already run the read.otu command
			if ((globaldata->gSparseMatrix == NULL) || (globaldata->gListVector == NULL)) {
				cout << "Before you use the get.oturep command, you first need to read in a distance matrix." << endl; 
				abort = true;
			}
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not found") { cout << "fasta is a required parameter for the get.oturep command." << endl; abort = true; }
			else if (fastafile == "not open") { abort = true; }	
		
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not found") { cout << "list is a required parameter for the get.oturep command." << endl; abort = true; }
			else if (listfile == "not open") { abort = true; }	

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
			else {
				//read in group map info.
				groupMap = new GroupMap(groupfile);
				groupMap->readMap();
			}
		
			if (abort == false) {
			
				if(globaldata->gSparseMatrix != NULL)	{	matrix = globaldata->gSparseMatrix;		}	
					
				//globaldata->gListVector bin 0 = first name read in distance matrix, globaldata->gListVector bin 1 = second name read in distance matrix
				if(globaldata->gListVector != NULL)		{
		
					vector<string> names;
					string binnames;
					//map names to rows in sparsematrix
					for (int i = 0; i < globaldata->gListVector->size(); i++) {
						names.clear();
						binnames = globaldata->gListVector->get(i);
	
						splitAtComma(binnames, names);
				
						for (int j = 0; j < names.size(); j++) {
							nameToIndex[names[j]] = i;
						}
	
					}
				}else { cout << "error, no listvector." << endl; }
				
				openInputFile(fastafile, in);
				fasta = new FastaMap();
			}
		
		}
	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetOTURepCommand class Function GetOTURepCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetOTURepCommand class function GetOTURepCommand. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
}

//**********************************************************************************************************************

void GetOTURepCommand::help(){
	try {
		cout << "The get.oturep command can only be executed after a successful read.dist command." << "\n";
		cout << "The get.oturep command parameters are list, fasta, name, group, line and label.  The fasta and list parameters are required, and you may not use line and label at the same time." << "\n";
		cout << "The line and label allow you to select what distance levels you would like a output files created for, and are separated by dashes." << "\n";
		cout << "The get.oturep command should be in the following format: get.oturep(fasta=yourFastaFile, list=yourListFile, name=yourNamesFile, group=yourGroupFile, line=yourLines, label=yourLabels)." << "\n";
		cout << "Example get.oturep(fasta=amazon.fasta, list=amazon.fn.list, group=amazon.groups, line=1-3-5, name=amazon.names)." << "\n";
		cout << "The default value for line and label are all lines in your inputfile." << "\n";
		cout << "The get.oturep command outputs a .fastarep file for each distance you specify, selecting one OTU representative for each bin." << "\n";
		cout << "If you provide a groupfile, then it also appends the names of the groups present in that bin." << "\n";
		cout << "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile)." << "\n" << "\n";
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetOTURepCommand class Function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetOTURepCommand class function help. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************

GetOTURepCommand::~GetOTURepCommand(){
	delete input;
	delete read;
	delete fasta;
	if (groupfile != "") {
		delete groupMap;
	}
}

//**********************************************************************************************************************

int GetOTURepCommand::execute(){
	try {
	
		if (abort == true) { return 0; }
		
		int count = 1;
		int error;
		
		//read fastafile
		fasta->readFastaFile(in);
		
		//set format to list so input can get listvector
//		globaldata->setFormat("list");
		
		//if user gave a namesfile then use it
		if (namesfile != "") {
			readNamesFile();
		}
		
		//read list file
		read = new ReadOTUFile(listfile);	
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
					cout << list->getLabel() << '\t' << count << endl;
					error = process(list);
					if (error == 1) { return 0; } //there is an error in hte input files, abort command
					
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
					userLines.erase(count);
			}
			
			if ((anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastList->getLabel()) != 1)) {
					cout << lastList->getLabel() << '\t' << count << endl;
					error = process(lastList);
					if (error == 1) { return 0; } //there is an error in hte input files, abort command
					
					processedLabels.insert(lastList->getLabel());
					userLabels.erase(lastList->getLabel());
			}
			
			if (count != 1) { delete lastList; }
			lastList = list;			
			
			list = input->getListVector();
			count++;
		}
		
		//output error messages about any remaining user labels
		bool needToRun = false;
		for (set<string>::iterator it = userLabels.begin(); it != userLabels.end(); it++) {  
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
			cout << lastList->getLabel() << '\t' << count << endl;
			error = process(lastList);
			if (error == 1) { return 0; } //there is an error in hte input files, abort command
		}
		delete lastList;
		
		delete matrix;
		globaldata->gSparseMatrix = NULL;
		delete list;
		globaldata->gListVector = NULL;

		return 0;
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetOTURepCommand class Function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetOTURepCommand class function execute. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}

}

//**********************************************************************************************************************
void GetOTURepCommand::readNamesFile() {
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
		cout << "Standard Error: " << e.what() << " has occurred in the GetOTURepCommand class Function readNamesFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetOTURepCommand class function readNamesFile. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}
//**********************************************************************************************************************
string GetOTURepCommand::FindRep(int bin, string& group, ListVector* thisList) {
	try{
		vector<string> names;
		map<string, float> sums;
		map<int, string> binMap; //subset of namesToIndex - just member of this bin
		string binnames;
		float min = 10000;
		string minName;
		map<string, string> groups;
		map<string, string>::iterator groupIt;
		
		binnames = thisList->get(bin);
	
		//parse names into vector
		splitAtComma(binnames, names);
		
		//if you have a groupfile
		if(groupfile != "") {
			//find the groups that are in this bin
			for (int i = 0; i < names.size(); i++) {
				string groupName = groupMap->getGroup(names[i]);
				if (groupName == "not found") {  
					cout << names[i] << " is missing from your group file. Please correct. " << endl;
					groupError = true;
				}else{
					groups[groupName] = groupName;
				}
			}
			
			//turn the groups into a string
			for(groupIt = groups.begin(); groupIt != groups.end(); groupIt++) { group += groupIt->first + "-"; }
			
			//rip off last dash
			group = group.substr(0, group.length()-1);
		}
		
		//if only 1 sequence in bin then that's the rep
		if (names.size() == 1) { return names[0]; }
		else {
			//fill binMap
			for (int i = 0; i < names.size(); i++) {
				for (map<string, int>::iterator it = nameToIndex.begin(); it != nameToIndex.end(); it++) {

					if (it->first == names[i]) {  
						binMap[it->second] = it->first;

						//initialize sums map
						sums[it->first] = 0.0;
						break;
					}
				}
			}
			
			//go through each cell in the sparsematrix
			for(MatData currentCell = matrix->begin(); currentCell != matrix->end(); currentCell++){
				//is this a distance between 2 members of this bin
				map<int, string>::iterator it = binMap.find(currentCell->row);
				map<int, string>::iterator it2 = binMap.find(currentCell->column);
				
				//sum the distance of the sequences in the bin to eachother
				if ((it != binMap.end()) && (it2 != binMap.end())) {
					//this is a cell that repesents the distance between to of this bins members
					sums[it->second] += currentCell->dist;
					sums[it2->second] += currentCell->dist;
				}
			}
			
			//smallest sum is the representative
			for (map<string, float>::iterator it4 = sums.begin(); it4 != sums.end(); it4++) {
				if (it4->second < min) {
					min = it4->second;
					minName = it4->first;
				}

			}
			
			return minName;
		}
	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetOTURepCommand class Function FindRep. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetOTURepCommand class function FindRep. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************
int GetOTURepCommand::process(ListVector* processList) {
	try{
				string nameRep, name, sequence;

				//create output file
				string outputFileName = getRootName(listfile) + processList->getLabel() + ".rep.fasta";
				openOutputFile(outputFileName, out);
				
				//for each bin in the list vector
				for (int i = 0; i < processList->size(); i++) {
					string groups;
					nameRep = FindRep(i, groups, processList);
					
					//print out name and sequence for that bin
					sequence = fasta->getSequence(nameRep);

					if (sequence != "not found") {
						if (groupfile == "") {
							nameRep = nameRep + "|" + toString(i+1);
							out << ">" << nameRep << endl;
							out << sequence << endl;
						}else {
							nameRep = nameRep + "|" + groups + "|" + toString(i+1);
							out << ">" << nameRep << endl;
							out << sequence << endl;
						}
					}else { 
						cout << nameRep << " is missing from your fasta or name file. Please correct. " << endl; 
						remove(outputFileName.c_str());
						return 1;
					}
				}
				
				out.close();
				return 0;
	
	}
	catch(exception& e) {
		cout << "Standard Error: " << e.what() << " has occurred in the GetOTURepCommand class Function process. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}
	catch(...) {
		cout << "An unknown error has occurred in the GetOTURepCommand class function process. Please contact Pat Schloss at pschloss@microbio.umass.edu." << "\n";
		exit(1);
	}	
}

//**********************************************************************************************************************





