/*
 *  getoturepcommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 4/6/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "getoturepcommand.h"
#include "readphylip.h"
#include "readcolumn.h"
#include "formatphylip.h"
#include "formatcolumn.h"



//********************************************************************************************************************
//sorts lowest to highest
inline bool compareName(repStruct left, repStruct right){
	return (left.name < right.name);	
}
//********************************************************************************************************************
//sorts lowest to highest
inline bool compareBin(repStruct left, repStruct right){
	return (left.bin < right.bin);	
}
//********************************************************************************************************************
//sorts lowest to highest
inline bool compareSize(repStruct left, repStruct right){
	return (left.size < right.size);	
}
//********************************************************************************************************************
//sorts lowest to highest
inline bool compareGroup(repStruct left, repStruct right){
	return (left.group < right.group);	
}
//**********************************************************************************************************************
GetOTURepCommand::GetOTURepCommand(string option)  {
	try{
		globaldata = GlobalData::getInstance();
		abort = false;
		allLines = 1;
		labels.clear();
				
		//allow user to run help
		if (option == "help") { 
			help(); abort = true;
		} else {
			//valid paramters for this command
			string Array[] =  {"fasta","list","label","name", "group", "sorted", "phylip","column","large","cutoff","precision","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string, string>::iterator it;
		
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
				
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("phylip");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["phylip"] = inputDir + it->second;		}
				}
				
				it = parameters.find("column");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["column"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
			}

			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not found") { m->mothurOut("fasta is a required parameter for the get.oturep command."); m->mothurOutEndLine(); abort = true; }
			else if (fastafile == "not open") { abort = true; }	
		
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not found") { m->mothurOut("list is a required parameter for the get.oturep command."); m->mothurOutEndLine(); abort = true; }
			else if (listfile == "not open") { abort = true; }	
			
			phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not found") { phylipfile = "";  }
			else if (phylipfile == "not open") { abort = true; }	
			else { distFile = phylipfile; format = "phylip";   }
			
			columnfile = validParameter.validFile(parameters, "column", true);
			if (columnfile == "not found") { columnfile = ""; }
			else if (columnfile == "not open") { abort = true; }	
			else { distFile = columnfile; format = "column";   }
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			
			if ((phylipfile == "") && (columnfile == "")) { m->mothurOut("When executing a get.oturep command you must enter a phylip or a column."); m->mothurOutEndLine(); abort = true; }
			else if ((phylipfile != "") && (columnfile != "")) { m->mothurOut("When executing a get.oturep command you must enter ONLY ONE of the following: phylip or column."); m->mothurOutEndLine(); abort = true; }
		
			if (columnfile != "") {  if (namefile == "") {  cout << "You need to provide a namefile if you are going to use the column format." << endl; abort = true; }  }

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; allLines = 1;  }
			else { 
				if(label != "all") {  splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { groupfile = ""; abort = true; }
			else if (groupfile == "not found") { groupfile = ""; }
			else {
				//read in group map info.
				groupMap = new GroupMap(groupfile);
				int error = groupMap->readMap();
				if (error == 1) { delete groupMap; abort = true; }
			}
			
			sorted = validParameter.validFile(parameters, "sorted", false);		if (sorted == "not found"){	sorted = "";	}
			if ((sorted != "") && (sorted != "name") && (sorted != "bin") && (sorted != "size") && (sorted != "group")) {
				m->mothurOut(sorted + " is not a valid option for the sorted parameter. The only options are: name, bin, size and group. I will not sort."); m->mothurOutEndLine();
				sorted = "";
			}
			
			if ((sorted == "group") && (groupfile == "")) {
				m->mothurOut("You must provide a groupfile to sort by group. I will not sort."); m->mothurOutEndLine();
				sorted = "";
			}
			
			string temp = validParameter.validFile(parameters, "large", false);		if (temp == "not found") {	temp = "F";	}
			large = isTrue(temp);
			
			temp = validParameter.validFile(parameters, "precision", false);			if (temp == "not found") { temp = "100"; }
			convert(temp, precision); 
			
			temp = validParameter.validFile(parameters, "cutoff", false);			if (temp == "not found") { temp = "10.0"; }
			convert(temp, cutoff); 
			cutoff += (5 / (precision * 10.0));
		}
	}
	catch(exception& e) {
		m->errorOut(e, "GetOTURepCommand", "GetOTURepCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

void GetOTURepCommand::help(){
	try {
		m->mothurOut("The get.oturep command parameters are phylip, column, list, fasta, name, group, large, cutoff, precision, sorted and label.  The fasta and list parameters are required, as well as phylip or column and name.\n");
		m->mothurOut("The label parameter allows you to select what distance levels you would like a output files created for, and is separated by dashes.\n");
		m->mothurOut("The phylip or column parameter is required, but only one may be used.  If you use a column file the name filename is required. \n");
		m->mothurOut("If you do not provide a cutoff value 10.00 is assumed. If you do not provide a precision value then 100 is assumed.\n");
		m->mothurOut("The get.oturep command should be in the following format: get.oturep(phylip=yourDistanceMatrix, fasta=yourFastaFile, list=yourListFile, name=yourNamesFile, group=yourGroupFile, label=yourLabels).\n");
		m->mothurOut("Example get.oturep(phylip=amazon.dist, fasta=amazon.fasta, list=amazon.fn.list, group=amazon.groups).\n");
		m->mothurOut("The default value for label is all labels in your inputfile.\n");
		m->mothurOut("The sorted parameter allows you to indicate you want the output sorted. You can sort by sequence name, bin number, bin size or group. The default is no sorting, but your options are name, number, size, or group.\n");
		m->mothurOut("The large parameter allows you to indicate that your distance matrix is too large to fit in RAM.  The default value is false.\n");
		m->mothurOut("The get.oturep command outputs a .fastarep and .rep.names file for each distance you specify, selecting one OTU representative for each bin.\n");
		m->mothurOut("If you provide a groupfile, then it also appends the names of the groups present in that bin.\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "GetOTURepCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

GetOTURepCommand::~GetOTURepCommand(){}

//**********************************************************************************************************************

int GetOTURepCommand::execute(){
	try {
	
		if (abort == true) { return 0; }
		int error;
		
		if (!large) {
			//read distance files
			if (format == "column") { readMatrix = new ReadColumnMatrix(distFile); }	
			else if (format == "phylip") { readMatrix = new ReadPhylipMatrix(distFile); }
			else { m->mothurOut("File format error."); m->mothurOutEndLine(); return 0;  }
			
			readMatrix->setCutoff(cutoff);
	
			if(namefile != ""){	
				nameMap = new NameAssignment(namefile);
				nameMap->readMap();
			}else{	nameMap = NULL;		}
			
			readMatrix->read(nameMap);
			
			if (m->control_pressed) { delete readMatrix; delete groupMap; return 0; }

			//get matrix
			if (globaldata->gListVector != NULL) {  delete globaldata->gListVector;  }
			globaldata->gListVector = readMatrix->getListVector();

			SparseMatrix* matrix = readMatrix->getMatrix();
			
			// Create a data structure to quickly access the distance information.
			// It consists of a vector of distance maps, where each map contains
			// all distances of a certain sequence. Vector and maps are accessed
			// via the index of a sequence in the distance matrix
			seqVec = vector<SeqMap>(globaldata->gListVector->size()); 
			for (MatData currentCell = matrix->begin(); currentCell != matrix->end(); currentCell++) {
				if (m->control_pressed) { delete readMatrix; delete groupMap; return 0; }
				seqVec[currentCell->row][currentCell->column] = currentCell->dist;
			}
			
			delete matrix;
			delete readMatrix;
			
			if (m->control_pressed) {  delete groupMap; return 0; }
		}else {
			//process file and set up indexes
			if (format == "column") { formatMatrix = new FormatColumnMatrix(distFile); }	
			else if (format == "phylip") { formatMatrix = new FormatPhylipMatrix(distFile); }
			else { m->mothurOut("File format error."); m->mothurOutEndLine(); return 0;  }
			
			formatMatrix->setCutoff(cutoff);
	
			if(namefile != ""){	
				nameMap = new NameAssignment(namefile);
				nameMap->readMap();
			}else{	nameMap = NULL;		}
			
			formatMatrix->read(nameMap);
			
			if (m->control_pressed) { delete formatMatrix; delete groupMap; return 0; }

			//get matrix
			if (globaldata->gListVector != NULL) {  delete globaldata->gListVector;  }
			globaldata->gListVector = formatMatrix->getListVector();
			
			distFile = formatMatrix->getFormattedFileName();
			
			//positions in file where the distances for each sequence begin
			//rowPositions[1] = position in file where distance related to sequence 1 start.
			rowPositions = formatMatrix->getRowPositions();
			
			delete formatMatrix;
			
			//openfile for getMap to use
			openInputFile(distFile, inRow);
			
			if (m->control_pressed) { inRow.close(); remove(distFile.c_str()); delete groupMap; return 0; }
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
		} else { m->mothurOut("error, no listvector."); m->mothurOutEndLine(); }
		
		fasta = new FastaMap();
		
		//read fastafile
		fasta->readFastaFile(fastafile);
		
		if (m->control_pressed) { 
			if (large) {  inRow.close(); remove(distFile.c_str());  }
			delete fasta; return 0; 
		}
				
		//if user gave a namesfile then use it
		if (namefile != "") {	readNamesFile();	}
		
		//set format to list so input can get listvector
		globaldata->setFormat("list");

		//read list file
		read = new ReadOTUFile(listfile);
		read->read(&*globaldata); 
		
		input = globaldata->ginput;
		list = globaldata->gListVector;
		string lastLabel = list->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		if (m->control_pressed) { 
			if (large) {  inRow.close(); remove(distFile.c_str());  }
			delete fasta; delete read; delete input; delete list; globaldata->gListVector = NULL; return 0; 
		}
	
		while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (allLines == 1 || labels.count(list->getLabel()) == 1){
					m->mothurOut(list->getLabel() + "\t" + toString(list->size())); m->mothurOutEndLine();
					error = process(list);
					if (error == 1) { return 0; } //there is an error in hte input files, abort command
					
					if (m->control_pressed) { 
						if (large) {  inRow.close(); remove(distFile.c_str());  }
						if (groupfile != "") {	delete groupMap;  globaldata->gGroupmap = NULL;	}
						for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  }
						delete fasta; delete read; delete input; delete list; globaldata->gListVector = NULL; return 0; 
					}
					
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
			}
			
			if ((anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = list->getLabel();
					
					delete list;
					list = input->getListVector(lastLabel);
					m->mothurOut(list->getLabel() + "\t" + toString(list->size())); m->mothurOutEndLine();
					error = process(list);
					if (error == 1) { return 0; } //there is an error in hte input files, abort command
					
					if (m->control_pressed) { 
						if (large) {  inRow.close(); remove(distFile.c_str());  }
						if (groupfile != "") {	delete groupMap;  globaldata->gGroupmap = NULL;	}
						for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  }
						delete fasta; delete read; delete input; delete list; globaldata->gListVector = NULL; return 0; 
					}
					
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
					
					//restore real lastlabel to save below
					list->setLabel(saveLabel);
			}
			
			lastLabel = list->getLabel();
			
			delete list;
			list = input->getListVector();
		}
		
		//output error messages about any remaining user labels
		bool needToRun = false;
		for (set<string>::iterator it = userLabels.begin(); it != userLabels.end(); it++) {  
			m->mothurOut("Your file does not include the label " + *it); 
			if (processedLabels.count(list->getLabel()) != 1) {
				m->mothurOut(". I will use " + lastLabel + "."); m->mothurOutEndLine();
				needToRun = true;
			}else {
				m->mothurOut(". Please refer to " + lastLabel + "."); m->mothurOutEndLine();
			}
		}
		
		//run last label if you need to
		if (needToRun == true)  {
			if (list != NULL) {	delete list;	}
			list = input->getListVector(lastLabel);
			m->mothurOut(list->getLabel() + "\t" + toString(list->size())); m->mothurOutEndLine();
			error = process(list);
			delete list;
			if (error == 1) { return 0; } //there is an error in hte input files, abort command
			
			if (m->control_pressed) { 
					if (large) {  inRow.close(); remove(distFile.c_str());  }
					if (groupfile != "") {	delete groupMap;  globaldata->gGroupmap = NULL;	}
					for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  }
					delete fasta; delete read; delete input; delete list; globaldata->gListVector = NULL; return 0; 
			}
		}
		
		//close and remove formatted matrix file
		if (large) {
			inRow.close();
			remove(distFile.c_str());
		}
		
		globaldata->gListVector = NULL;
		delete input;  globaldata->ginput = NULL;
		delete read;
		delete fasta;
		if (groupfile != "") {
			delete groupMap;  globaldata->gGroupmap = NULL;
		}
		
		if (m->control_pressed) {  return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "GetOTURepCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
void GetOTURepCommand::readNamesFile() {
	try {
		vector<string> dupNames;
		openInputFile(namefile, inNames);
		
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
		m->errorOut(e, "GetOTURepCommand", "readNamesFile");
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
					m->mothurOut(names[i] + " is missing from your group file. Please correct. "); m->mothurOutEndLine();
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
		}else{ group = ""; }

		// if only 1 sequence in bin or processing the "unique" label, then 
		// the first sequence of the OTU is the representative one
		if ((names.size() == 1) || (list->getLabel() == "unique")) {
			return names[0];
		}else{
			vector<int> seqIndex(names.size());
			vector<float> max_dist(names.size());
			vector<float> total_dist(names.size());

			//fill seqIndex and initialize sums
			for (size_t i = 0; i < names.size(); i++) {
				seqIndex[i] = nameToIndex[names[i]];
				max_dist[i] = 0.0;
				total_dist[i] = 0.0;
			}
			
			// loop through all entries in seqIndex
			SeqMap::iterator it;
			SeqMap currMap;
			for (size_t i=0; i < seqIndex.size(); i++) {
				if (m->control_pressed) {  return  "control"; }
			
				if (!large) {	currMap = seqVec[seqIndex[i]];  }
				else		{	currMap = getMap(seqIndex[i]);	}
				
				for (size_t j=0; j < seqIndex.size(); j++) {
					it = currMap.find(seqIndex[j]);		
					if (it != currMap.end()) {
						max_dist[i] = max(max_dist[i], it->second);
						max_dist[j] = max(max_dist[j], it->second);
						total_dist[i] += it->second;
						total_dist[j] += it->second;
					}else{ //if you can't find the distance make it the cutoff
						max_dist[i] = max(max_dist[i], cutoff);
						max_dist[j] = max(max_dist[j], cutoff);
						total_dist[i] += cutoff;
						total_dist[j] += cutoff;
					}
				}
			}
			
			// sequence with the smallest maximum distance is the representative
			//if tie occurs pick sequence with smallest average distance
			float min = 10000;
			int minIndex;
			for (size_t i=0; i < max_dist.size(); i++) {
				if (m->control_pressed) {  return  "control"; }
				if (max_dist[i] < min) {
					min = max_dist[i];
					minIndex = i;
				}else if (max_dist[i] == min) {
					float currentAverage = total_dist[minIndex] / (float) total_dist.size();
					float newAverage = total_dist[i] / (float) total_dist.size();
					
					if (newAverage < currentAverage) {
						min = max_dist[i];
						minIndex = i;
					}
				}
			}

			return(names[minIndex]);
		}
	}
	catch(exception& e) {
		m->errorOut(e, "GetOTURepCommand", "FindRep");
		exit(1);
	}
}

//**********************************************************************************************************************
int GetOTURepCommand::process(ListVector* processList) {
	try{
		string nameRep, name, sequence;

		//create output file
		if (outputDir == "") { outputDir += hasPath(listfile); }
		string outputFileName = outputDir + getRootName(getSimpleName(listfile)) + processList->getLabel() + ".rep.fasta";
		openOutputFile(outputFileName, out);
		vector<repStruct> reps;
		outputNames.push_back(outputFileName);
		
		ofstream newNamesOutput;
		string outputNamesFile = outputDir + getRootName(getSimpleName(listfile)) + processList->getLabel() + ".rep.names";
		openOutputFile(outputNamesFile, newNamesOutput);
		outputNames.push_back(outputNamesFile);
		
		//for each bin in the list vector
		for (int i = 0; i < processList->size(); i++) {
			string groups;
			int binsize;
			
			if (m->control_pressed) { out.close();  newNamesOutput.close(); return 0; }
			
			nameRep = findRep(i, groups, processList, binsize);
			
			if (m->control_pressed) { out.close();  newNamesOutput.close(); return 0; }
			
			//output to new names file
			newNamesOutput << nameRep << '\t' << processList->get(i) << endl;

			//print out name and sequence for that bin
			sequence = fasta->getSequence(nameRep);

			if (sequence != "not found") {
				if (sorted == "") { //print them out
					nameRep = nameRep + "|" + toString(i+1);
					nameRep = nameRep + "|" + toString(binsize);
					if (groupfile != "") {
						nameRep = nameRep + "|" + groups;
					}
					out << ">" << nameRep << endl;
					out << sequence << endl;
				}else { //save them
					repStruct newRep(nameRep, i+1, binsize, groups);
					reps.push_back(newRep);
				}
			}else { 
				m->mothurOut(nameRep + " is missing from your fasta or name file. Please correct. "); m->mothurOutEndLine(); 
				remove(outputFileName.c_str());
				remove(outputNamesFile.c_str());
				return 1;
			}
		}
		
		if (sorted != "") { //then sort them and print them
			if (sorted == "name")		{  sort(reps.begin(), reps.end(), compareName);		}
			else if (sorted == "bin")	{  sort(reps.begin(), reps.end(), compareBin);		}
			else if (sorted == "size")	{  sort(reps.begin(), reps.end(), compareSize);		}
			else if (sorted == "group")	{  sort(reps.begin(), reps.end(), compareGroup);	}
			
			//print them
			for (int i = 0; i < reps.size(); i++) {
				string sequence = fasta->getSequence(reps[i].name);
				string outputName = reps[i].name + "|" + toString(reps[i].bin);
				outputName = outputName + "|" + toString(reps[i].size);
				if (groupfile != "") {
					outputName = outputName + "|" + reps[i].group;
				}
				out << ">" << outputName << endl;
				out << sequence << endl;
			}
		}

		out.close();
		newNamesOutput.close();
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "GetOTURepCommand", "process");
		exit(1);
	}
}

//**********************************************************************************************************************
SeqMap GetOTURepCommand::getMap(int row) {
	try {
		SeqMap rowMap;
		
		//make sure this row exists in the file, it may not if the seq did not have any distances below the cutoff
		if (rowPositions[row] != -1){
			//go to row in file
			inRow.seekg(rowPositions[row]);
			
			int rowNum, numDists, colNum;
			float dist;
			
			inRow >> rowNum >> numDists;
			
			for(int i = 0; i < numDists; i++) {
				inRow >> colNum >> dist;
				rowMap[colNum] = dist;
				
			}
		}
		
		return rowMap;
	}
	catch(exception& e) {
		m->errorOut(e, "GetOTURepCommand", "getMap");
		exit(1);
	}
}
//**********************************************************************************************************************

