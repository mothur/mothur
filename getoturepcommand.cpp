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
		if (option == "help") { 
			help(); abort = true;
		} else {
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
				mothurOut("Before you use the get.oturep command, you first need to read in a distance matrix."); mothurOutEndLine(); 
				abort = true;
			}
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not found") { mothurOut("fasta is a required parameter for the get.oturep command."); mothurOutEndLine(); abort = true; }
			else if (fastafile == "not open") { abort = true; }	
		
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not found") { mothurOut("list is a required parameter for the get.oturep command."); mothurOutEndLine(); abort = true; }
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
			if ((line != "") && (label != "")) { mothurOut("You cannot use both the line and label parameters at the same time. "); mothurOutEndLine(); abort = true; }
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
			
				if(globaldata->gSparseMatrix != NULL)	{
					matrix = globaldata->gSparseMatrix;
				}

				//globaldata->gListVector bin 0 = first name read in distance matrix, globaldata->gListVector bin 1 = second name read in distance matrix
				if (globaldata->gListVector != NULL) {
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
				} else { mothurOut("error, no listvector."); mothurOutEndLine(); }

				// Create a data structure to quickly access the distance information.
				// It consists of a vector of distance maps, where each map contains
				// all distances of a certain sequence. Vector and maps are accessed
				// via the index of a sequence in the distance matrix
				seqVec = vector<SeqMap>(globaldata->gListVector->size()); 
				for (MatData currentCell = matrix->begin(); currentCell != matrix->end(); currentCell++) {
					seqVec[currentCell->row][currentCell->column] = currentCell->dist;
				}

				fasta = new FastaMap();
			}
		}
	}
	catch(exception& e) {
		errorOut(e, "GetOTURepCommand", "GetOTURepCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void GetOTURepCommand::help(){
	try {
		mothurOut("The get.oturep command can only be executed after a successful read.dist command.\n");
		mothurOut("The get.oturep command parameters are list, fasta, name, group, line and label.  The fasta and list parameters are required, and you may not use line and label at the same time.\n");
		mothurOut("The line and label allow you to select what distance levels you would like a output files created for, and are separated by dashes.\n");
		mothurOut("The get.oturep command should be in the following format: get.oturep(fasta=yourFastaFile, list=yourListFile, name=yourNamesFile, group=yourGroupFile, line=yourLines, label=yourLabels).\n");
		mothurOut("Example get.oturep(fasta=amazon.fasta, list=amazon.fn.list, group=amazon.groups, line=1-3-5, name=amazon.names).\n");
		mothurOut("The default value for line and label are all lines in your inputfile.\n");
		mothurOut("The get.oturep command outputs a .fastarep file for each distance you specify, selecting one OTU representative for each bin.\n");
		mothurOut("If you provide a groupfile, then it also appends the names of the groups present in that bin.\n");
		mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n\n");
	}
	catch(exception& e) {
		errorOut(e, "GetOTURepCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

GetOTURepCommand::~GetOTURepCommand(){
	if (abort == false) {
		delete input;  globaldata->ginput = NULL;
		delete read;
		delete fasta;
		if (groupfile != "") {
			delete groupMap;  globaldata->gGroupmap = NULL;
		}
	}
}

//**********************************************************************************************************************

int GetOTURepCommand::execute(){
	try {
	
		if (abort == true) { return 0; }

		int count = 1;
		int error;
		
		//read fastafile
		fasta->readFastaFile(fastafile);
		
		in.close();
		
		//set format to list so input can get listvector
		globaldata->setFormat("list");
		
		//if user gave a namesfile then use it
		if (namesfile != "") {
			readNamesFile();
		}
		
		//read list file
		read = new ReadOTUFile(listfile);
		read->read(&*globaldata); 
		
		input = globaldata->ginput;
		list = globaldata->gListVector;
		string lastLabel = list->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		set<int> userLines = lines;

		
		while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0) || (userLines.size() != 0))) {
			
			if (allLines == 1 || lines.count(count) == 1 || labels.count(list->getLabel()) == 1){
					mothurOut(list->getLabel() + "\t" + toString(list->size())); mothurOutEndLine();
					error = process(list);
					if (error == 1) { return 0; } //there is an error in hte input files, abort command
					
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
					userLines.erase(count);
			}
			
			if ((anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					delete list;
					list = input->getListVector(lastLabel);
					mothurOut(list->getLabel() + "\t" + toString(list->size())); mothurOutEndLine();
					error = process(list);
					if (error == 1) { return 0; } //there is an error in hte input files, abort command
					
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
			}
			
			lastLabel = list->getLabel();
			
			delete list;
			list = input->getListVector();
			count++;
		}
		
		//output error messages about any remaining user labels
		bool needToRun = false;
		for (set<string>::iterator it = userLabels.begin(); it != userLabels.end(); it++) {  
			mothurOut("Your file does not include the label " + *it); 
			if (processedLabels.count(list->getLabel()) != 1) {
				mothurOut(". I will use " + lastLabel + "."); mothurOutEndLine();
				needToRun = true;
			}else {
				mothurOut(". Please refer to " + lastLabel + "."); mothurOutEndLine();
			}
		}
		
		//run last line if you need to
		if (needToRun == true)  {
			if (list != NULL) {	delete list;	}
			list = input->getListVector(lastLabel);
			mothurOut(list->getLabel() + "\t" + toString(list->size())); mothurOutEndLine();
			error = process(list);
			delete list;
			if (error == 1) { return 0; } //there is an error in hte input files, abort command
		}
		
		delete matrix;
		globaldata->gSparseMatrix = NULL;
		delete list;
		globaldata->gListVector = NULL;

		return 0;
	}
	catch(exception& e) {
		errorOut(e, "GetOTURepCommand", "execute");
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
		errorOut(e, "GetOTURepCommand", "readNamesFile");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetOTURepCommand::findRep(int bin, string& group, ListVector* thisList, int& binsize) {
	try{
		vector<string> names;
		map<string, string> groups;
		map<string, string>::iterator groupIt;

		//parse names into vector
		string binnames = thisList->get(bin);
		splitAtComma(binnames, names);
		binsize = names.size();

		//if you have a groupfile
		if (groupfile != "") {
			//find the groups that are in this bin
			for (size_t i = 0; i < names.size(); i++) {
				string groupName = groupMap->getGroup(names[i]);
				if (groupName == "not found") {  
					mothurOut(names[i] + " is missing from your group file. Please correct. "); mothurOutEndLine();
					groupError = true;
				} else {
					groups[groupName] = groupName;
				}
			}
			
			//turn the groups into a string
			for (groupIt = groups.begin(); groupIt != groups.end(); groupIt++) {
				group += groupIt->first + "-";
			}
			//rip off last dash
			group = group.substr(0, group.length()-1);
		}

		// if only 1 sequence in bin or processing the "unique" line, then 
		// the first sequence of the OTU is the representative one
		if ((names.size() == 1) || (list->getLabel() == "unique")) {
			return names[0];
		} else {
			vector<int> seqIndex(names.size());
			vector<float> max_dist(names.size());

			//fill seqIndex and initialize sums
			for (size_t i = 0; i < names.size(); i++) {
				seqIndex[i] = nameToIndex[names[i]];
				max_dist[i] = 0.0;
			}
			
			// loop through all entries in seqIndex
			SeqMap::iterator it;
			SeqMap currMap;
			for (size_t i=0; i < seqIndex.size(); i++) {
				currMap = seqVec[seqIndex[i]];
				for (size_t j=0; j < seqIndex.size(); j++) {
					it = currMap.find(seqIndex[j]);		
					if (it != currMap.end()) {
						max_dist[i] = max(max_dist[i], it->second);
						max_dist[j] = max(max_dist[j], it->second);
					}
				}
			}
			
			// sequence with the smallest maximum distance is the representative
			float min = 10000;
			int minIndex;
			for (size_t i=0; i < max_dist.size(); i++) {
				if (max_dist[i] < min) {
					min = max_dist[i];
					minIndex = i;
				}
			}

			return(names[minIndex]);
		}
	}
	catch(exception& e) {
		errorOut(e, "GetOTURepCommand", "FindRep");
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
			int binsize;
			nameRep = findRep(i, groups, processList, binsize);

			//print out name and sequence for that bin
			sequence = fasta->getSequence(nameRep);

			if (sequence != "not found") {
				nameRep = nameRep + "|" + toString(i+1);
				nameRep = nameRep + "|" + toString(binsize);
				if (groupfile != "") {
					nameRep = nameRep + "|" + groups;
				}
				out << ">" << nameRep << endl;
				out << sequence << endl;
			} else { 
				mothurOut(nameRep + " is missing from your fasta or name file. Please correct. "); mothurOutEndLine(); 
				remove(outputFileName.c_str());
				return 1;
			}
		}

		out.close();
		return 0;

	}
	catch(exception& e) {
		errorOut(e, "GetOTURepCommand", "process");
		exit(1);
	}
}

//**********************************************************************************************************************
