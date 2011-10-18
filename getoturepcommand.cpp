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
#include "sharedutilities.h"


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
vector<string> GetOTURepCommand::setParameters(){	
	try {
		CommandParameter plist("list", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(plist);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(pfasta);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "none", "none",false,false); parameters.push_back(pgroup);
		CommandParameter pphylip("phylip", "InputTypes", "", "", "PhylipColumn", "PhylipColumn", "none",false,false); parameters.push_back(pphylip);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "ColumnName",false,false); parameters.push_back(pname);
		CommandParameter pcolumn("column", "InputTypes", "", "", "PhylipColumn", "PhylipColumn", "ColumnName",false,false); parameters.push_back(pcolumn);
		CommandParameter plabel("label", "String", "", "", "", "", "",false,false); parameters.push_back(plabel);
		CommandParameter pgroups("groups", "String", "", "", "", "", "",false,false); parameters.push_back(pgroups);
		CommandParameter pcutoff("cutoff", "Number", "", "10", "", "", "",false,false); parameters.push_back(pcutoff);
		CommandParameter pprecision("precision", "Number", "", "100", "", "", "",false,false); parameters.push_back(pprecision);
		CommandParameter pweighted("weighted", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(pweighted);
		CommandParameter psorted("sorted", "Multiple", "none-name-bin-size-group", "none", "", "", "",false,false); parameters.push_back(psorted);
		CommandParameter plarge("large", "Boolean", "", "F", "", "", "",false,false); parameters.push_back(plarge);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetOTURepCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetOTURepCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The get.oturep command parameters are phylip, column, list, fasta, name, group, large, weighted, cutoff, precision, groups, sorted and label.  The fasta and list parameters are required, as well as phylip or column and name, unless you have valid current files.\n";
		helpString += "The label parameter allows you to select what distance levels you would like a output files created for, and is separated by dashes.\n";
		helpString += "The phylip or column parameter is required, but only one may be used.  If you use a column file the name filename is required. \n";
		helpString += "If you do not provide a cutoff value 10.00 is assumed. If you do not provide a precision value then 100 is assumed.\n";
		helpString += "The get.oturep command should be in the following format: get.oturep(phylip=yourDistanceMatrix, fasta=yourFastaFile, list=yourListFile, name=yourNamesFile, group=yourGroupFile, label=yourLabels).\n";
		helpString += "Example get.oturep(phylip=amazon.dist, fasta=amazon.fasta, list=amazon.fn.list, group=amazon.groups).\n";
		helpString += "The default value for label is all labels in your inputfile.\n";
		helpString += "The sorted parameter allows you to indicate you want the output sorted. You can sort by sequence name, bin number, bin size or group. The default is no sorting, but your options are name, number, size, or group.\n";
		helpString += "The large parameter allows you to indicate that your distance matrix is too large to fit in RAM.  The default value is false.\n";
		helpString += "The weighted parameter allows you to indicate that want to find the weighted representative. You must provide a namesfile to set weighted to true.  The default value is false.\n";
		helpString += "The representative is found by selecting the sequence that has the smallest total distance to all other sequences in the OTU. If a tie occurs the smallest average distance is used.\n";
		helpString += "For weighted = false, mothur assumes the distance file contains only unique sequences, the list file may contain all sequences, but only the uniques are considered to become the representative. If your distance file contains all the sequences it would become weighted=true.\n";
		helpString += "For weighted = true, mothur assumes the distance file contains only unique sequences, the list file must contain all sequences, all sequences are considered to become the representative, but unique name will be used in the output for consistency.\n";
		helpString += "If your distance file contains all the sequence and you do not provide a name file, the weighted representative will be given, unless your listfile is unique. If you provide a namefile, then you can select weighted or unweighted.\n";
		helpString += "The group parameter allows you provide a group file.\n";
		helpString += "The groups parameter allows you to indicate that you want representative sequences for each group specified for each OTU, group name should be separated by dashes. ex. groups=A-B-C.\n";
		helpString += "The get.oturep command outputs a .fastarep and .rep.names file for each distance you specify, selecting one OTU representative for each bin.\n";
		helpString += "If you provide a groupfile, then it also appends the names of the groups present in that bin.\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "GetOTURepCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
GetOTURepCommand::GetOTURepCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["name"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "GetOTURepCommand", "GetOTURepCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
GetOTURepCommand::GetOTURepCommand(string option)  {
	try{
		abort = false; calledHelp = false;   
		allLines = 1;
				
		//allow user to run help
		if (option == "help") { 
			help(); abort = true; calledHelp = true;
		}else if(option == "citation") { citation(); abort = true; calledHelp = true;
		} else {
			vector<string> myArray = setParameters();
			
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
			outputTypes["name"] = tempOutNames;
			
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
				
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
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

			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not found") { 				
				fastafile = m->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (fastafile == "not open") { abort = true; }	
			else { m->setFastaFile(fastafile); }
		
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not found") { 			
				listfile = m->getListFile(); 
				if (listfile != "") { m->mothurOut("Using " + listfile + " as input file for the list parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current list file and the list parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (listfile == "not open") { abort = true; }	
			else { m->setListFile(listfile); }
			
			phylipfile = validParameter.validFile(parameters, "phylip", true);
			if (phylipfile == "not found") { phylipfile = "";  }
			else if (phylipfile == "not open") { abort = true; }	
			else { distFile = phylipfile; format = "phylip"; m->setPhylipFile(phylipfile);   }
			
			columnfile = validParameter.validFile(parameters, "column", true);
			if (columnfile == "not found") { columnfile = ""; }
			else if (columnfile == "not open") { abort = true; }	
			else { distFile = columnfile; format = "column";  m->setColumnFile(columnfile); }
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			else { m->setNameFile(namefile); }
			
			if ((phylipfile == "") && (columnfile == "")) { //is there are current file available for either of these?
				//give priority to column, then phylip
				columnfile = m->getColumnFile(); 
				if (columnfile != "") {  distFile = columnfile; format = "column"; m->mothurOut("Using " + columnfile + " as input file for the column parameter."); m->mothurOutEndLine(); }
				else { 
					phylipfile = m->getPhylipFile(); 
					if (phylipfile != "") {  distFile = phylipfile; format = "phylip"; m->mothurOut("Using " + phylipfile + " as input file for the phylip parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("No valid current files. You must provide a phylip or column file before you can use the get.oturep command."); m->mothurOutEndLine(); 
						abort = true;
					}
				}
			}else if ((phylipfile != "") && (columnfile != "")) { m->mothurOut("When executing a get.oturep command you must enter ONLY ONE of the following: phylip or column."); m->mothurOutEndLine(); abort = true; }
		
			if (columnfile != "") {  
				if (namefile == "") {  
					namefile = m->getNameFile(); 
					if (namefile != "") {  m->mothurOut("Using " + namefile + " as input file for the name parameter."); m->mothurOutEndLine(); }
					else { 
						m->mothurOut("You need to provide a namefile if you are going to use the column format."); m->mothurOutEndLine(); 
						abort = true; 
					}	
				} 
			}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = ""; allLines = 1;  }
			else { 
				if(label != "all") {  m->splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { groupfile = ""; abort = true; }
			else if (groupfile == "not found") { groupfile = ""; }
			else { m->setGroupFile(groupfile); }
			
			sorted = validParameter.validFile(parameters, "sorted", false);		if (sorted == "not found"){	sorted = "";	}
			if (sorted == "none") { sorted=""; }
			if ((sorted != "") && (sorted != "name") && (sorted != "bin") && (sorted != "size") && (sorted != "group")) {
				m->mothurOut(sorted + " is not a valid option for the sorted parameter. The only options are: name, bin, size and group. I will not sort."); m->mothurOutEndLine();
				sorted = "";
			}
			
			if ((sorted == "group") && (groupfile == "")) {
				m->mothurOut("You must provide a groupfile to sort by group. I will not sort."); m->mothurOutEndLine();
				sorted = "";
			}
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { groups = ""; }
			else { 
				if (groupfile == "") {
					m->mothurOut("You must provide a groupfile to use groups."); m->mothurOutEndLine();
					abort = true;
				}else { 
					m->splitAtDash(groups, Groups);
				}
			}
			m->setGroups(Groups);
			
			string temp = validParameter.validFile(parameters, "large", false);		if (temp == "not found") {	temp = "F";	}
			large = m->isTrue(temp);
			
			temp = validParameter.validFile(parameters, "weighted", false);		if (temp == "not found") {	 temp = "f"; 	}
			weighted = m->isTrue(temp);
			
			if ((weighted) && (namefile == "")) { m->mothurOut("You cannot set weighted to true unless you provide a namesfile."); m->mothurOutEndLine(); abort = true; }
			
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

int GetOTURepCommand::execute(){
	try {
	
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		int error;
		list = NULL;
		
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
			
			if (m->control_pressed) { delete readMatrix; return 0; }

			list = readMatrix->getListVector();

			SparseMatrix* matrix = readMatrix->getMatrix();
			
			// Create a data structure to quickly access the distance information.
			// It consists of a vector of distance maps, where each map contains
			// all distances of a certain sequence. Vector and maps are accessed
			// via the index of a sequence in the distance matrix
			seqVec = vector<SeqMap>(list->size()); 
			for (MatData currentCell = matrix->begin(); currentCell != matrix->end(); currentCell++) {
				if (m->control_pressed) { delete readMatrix; return 0; }
				seqVec[currentCell->row][currentCell->column] = currentCell->dist;
			}
			//add dummy map for unweighted calc
			SeqMap dummy;
			seqVec.push_back(dummy);
			
			delete matrix;
			delete readMatrix;
			delete nameMap;
			
			if (m->control_pressed) { return 0; }
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
			
			if (m->control_pressed) { delete formatMatrix;  return 0; }

			list = formatMatrix->getListVector();
			
			distFile = formatMatrix->getFormattedFileName();
			
			//positions in file where the distances for each sequence begin
			//rowPositions[1] = position in file where distance related to sequence 1 start.
			rowPositions = formatMatrix->getRowPositions();
			rowPositions.push_back(-1); //dummy row for unweighted calc
			
			delete formatMatrix;
			delete nameMap;
			
			//openfile for getMap to use
			m->openInputFile(distFile, inRow);
			
			if (m->control_pressed) { inRow.close(); m->mothurRemove(distFile); return 0; }
		}
		
		
		//list bin 0 = first name read in distance matrix, list bin 1 = second name read in distance matrix
		if (list != NULL) {
			vector<string> names;
			string binnames;
			//map names to rows in sparsematrix
			for (int i = 0; i < list->size(); i++) {
				names.clear();
				binnames = list->get(i);
				
				m->splitAtComma(binnames, names);
				
				for (int j = 0; j < names.size(); j++) {
					nameToIndex[names[j]] = i;
				}
			}
		} else { m->mothurOut("error, no listvector."); m->mothurOutEndLine(); }
		
				
		if (m->control_pressed) { 
			if (large) {  inRow.close(); m->mothurRemove(distFile);  }
			return 0; 
		}
		
		if (groupfile != "") {
			//read in group map info.
			groupMap = new GroupMap(groupfile);
			int error = groupMap->readMap();
			if (error == 1) { delete groupMap; m->mothurOut("Error reading your groupfile. Proceeding without groupfile."); m->mothurOutEndLine(); groupfile = "";  }
			
			if (Groups.size() != 0) {
				SharedUtil* util = new SharedUtil();
				vector<string> gNamesOfGroups = groupMap->getNamesOfGroups();
				util->setGroups(Groups, gNamesOfGroups, "getoturep");
				groupMap->setNamesOfGroups(gNamesOfGroups);
				delete util;
			}
		}
		
		//done with listvector from matrix
		if (list != NULL) { delete list; }
		
		input = new InputData(listfile, "list");
		list = input->getListVector();
		string lastLabel = list->getLabel();

		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		if (m->control_pressed) { 
			if (large) {  inRow.close(); m->mothurRemove(distFile);  }
			delete input; delete list; return 0; 
		}
		
		if ((!weighted) && (namefile != "")) { readNamesFile(weighted); }
		
		while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
			if (allLines == 1 || labels.count(list->getLabel()) == 1){
					m->mothurOut(list->getLabel() + "\t" + toString(list->size())); m->mothurOutEndLine();
					error = process(list);
					if (error == 1) { return 0; } //there is an error in hte input files, abort command
					
					if (m->control_pressed) { 
						if (large) {  inRow.close(); m->mothurRemove(distFile);  }
						for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } outputTypes.clear();
						delete input; delete list; return 0; 
					}
					
					processedLabels.insert(list->getLabel());
					userLabels.erase(list->getLabel());
			}
			
			if ((m->anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
					string saveLabel = list->getLabel();
					
					delete list;
					list = input->getListVector(lastLabel);
					m->mothurOut(list->getLabel() + "\t" + toString(list->size())); m->mothurOutEndLine();
					error = process(list);
					if (error == 1) { return 0; } //there is an error in hte input files, abort command
					
					if (m->control_pressed) { 
						if (large) {  inRow.close(); m->mothurRemove(distFile);  }
						for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } outputTypes.clear();
						delete input; delete list; return 0; 
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
			m->mothurOut("Your file does not include the label " + (*it)); 
			if (processedLabels.count(lastLabel) != 1) {
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
					if (large) {  inRow.close(); m->mothurRemove(distFile);  }
					for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } outputTypes.clear();
					delete input; delete list; return 0; 
			}
		}
		
		//close and remove formatted matrix file
		if (large) {
			inRow.close();
			m->mothurRemove(distFile);
		}
		
		delete input;  
		
		if (!weighted) { nameFileMap.clear(); }
		
		//read fastafile
		fasta = new FastaMap();
		fasta->readFastaFile(fastafile);
		
		//if user gave a namesfile then use it
		if (namefile != "") {	readNamesFile();	}
		
		//output create and output the .rep.fasta files
		map<string, string>::iterator itNameFile;
		for (itNameFile = outputNameFiles.begin(); itNameFile != outputNameFiles.end(); itNameFile++) {
			processNames(itNameFile->first, itNameFile->second);
		}
		
		delete fasta;
		if (groupfile != "") { delete groupMap;  }
		
		if (m->control_pressed) {  return 0; }
		
		//set fasta file as new current fastafile - use first one??
		string current = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
		}
		
		itTypes = outputTypes.find("name");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setNameFile(current); }
		}
		
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
		ifstream in;
		vector<string> dupNames;
		m->openInputFile(namefile, in);
		
		string name, names, sequence;
	
		while(!in.eof()){
			in >> name;			//read from first column  A
			in >> names;		//read from second column  A,B,C,D
			
			dupNames.clear();
			
			//parse names into vector
			m->splitAtComma(names, dupNames);
			
			//store names in fasta map
			sequence = fasta->getSequence(name);
			for (int i = 0; i < dupNames.size(); i++) {
				fasta->push_back(dupNames[i], sequence);
			}
		
			m->gobble(in);
		}
		in.close();

	}
	catch(exception& e) {
		m->errorOut(e, "GetOTURepCommand", "readNamesFile");
		exit(1);
	}
}
//**********************************************************************************************************************
//read names file to find the weighted rep for each bin
void GetOTURepCommand::readNamesFile(bool w) {
	try {
		ifstream in;
		vector<string> dupNames;
		m->openInputFile(namefile, in);
		
		string name, names, sequence;
		
		while(!in.eof()){
			in >> name;	m->gobble(in);		//read from first column  A
			in >> names;							//read from second column  A,B,C,D
			
			dupNames.clear();
			
			//parse names into vector
			m->splitAtComma(names, dupNames);
			
			for (int i = 0; i < dupNames.size(); i++) {
				nameFileMap[dupNames[i]] = name;
			}
			
			m->gobble(in);
		}
		in.close();
		
	}
	catch(exception& e) {
		m->errorOut(e, "GetOTURepCommand", "readNamesFile");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetOTURepCommand::findRep(vector<string> names) {
	try{
		// if only 1 sequence in bin or processing the "unique" label, then 
		// the first sequence of the OTU is the representative one
		if ((names.size() == 1)) {
			return names[0];
		}else{
			vector<int> seqIndex(names.size());
			vector<float> max_dist(names.size());
			vector<float> total_dist(names.size());
			map<string, string>::iterator itNameFile;
			map<string, int>::iterator itNameIndex;

			//fill seqIndex and initialize sums
			for (size_t i = 0; i < names.size(); i++) {
				if (weighted) {
					seqIndex[i] = nameToIndex[names[i]];
				}else { 
					if (namefile == "") {
						itNameIndex = nameToIndex.find(names[i]);
						
						if (itNameIndex == nameToIndex.end()) { // you are not in the distance file and no namesfile, then assume you are not unique
							if (large) {  seqIndex[i] = (rowPositions.size()-1); }
							else {  seqIndex[i] = (seqVec.size()-1); }
						}else {
							seqIndex[i] = itNameIndex->second;
						}
						
					}else {
						itNameFile = nameFileMap.find(names[i]);
						
						if (itNameFile == nameFileMap.end()) {
							m->mothurOut("[ERROR]: " + names[i] + " is not in your namefile, please correct."); m->mothurOutEndLine(); m->control_pressed = true; 
						}else{
							string name1 = itNameFile->first;
							string name2 = itNameFile->second;
							
							if (name1 == name2) { //then you are unique so add your real dists
								seqIndex[i] = nameToIndex[names[i]];
							}else { //add dummy
								if (large) {  seqIndex[i] = (rowPositions.size()-1); }
								else {  seqIndex[i] = (seqVec.size()-1); }
							}
						}
					}
				}
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
		string name, sequence;
		string nameRep;

		//create output file
		if (outputDir == "") { outputDir += m->hasPath(listfile); }
				
		ofstream newNamesOutput;
		string outputNamesFile;
		map<string, ofstream*> filehandles;
		
		if (Groups.size() == 0) { //you don't want to use groups
			outputNamesFile  = outputDir + m->getRootName(m->getSimpleName(listfile)) + processList->getLabel() + ".rep.names";
			m->openOutputFile(outputNamesFile, newNamesOutput);
			outputNames.push_back(outputNamesFile); outputTypes["name"].push_back(outputNamesFile); 
			outputNameFiles[outputNamesFile] = processList->getLabel();
		}else{ //you want to use groups
			ofstream* temp;
			for (int i=0; i<Groups.size(); i++) {
				temp = new ofstream;
				filehandles[Groups[i]] = temp;
				outputNamesFile = outputDir + m->getRootName(m->getSimpleName(listfile)) + processList->getLabel() + "." + Groups[i] + ".rep.names";
				
				m->openOutputFile(outputNamesFile, *(temp));
				outputNames.push_back(outputNamesFile); outputTypes["name"].push_back(outputNamesFile);
				outputNameFiles[outputNamesFile] = processList->getLabel() + "." + Groups[i];
			}
		}
		
		//for each bin in the list vector
		for (int i = 0; i < processList->size(); i++) {
			if (m->control_pressed) { 
				out.close();  
				if (Groups.size() == 0) { //you don't want to use groups
					newNamesOutput.close();
				}else{
					for (int j=0; j<Groups.size(); j++) {
						(*(filehandles[Groups[j]])).close();
						delete filehandles[Groups[j]];
					}
				}
				return 0; 
			}
			
			string temp = processList->get(i);
			vector<string> namesInBin;
			m->splitAtComma(temp, namesInBin);
			
			if (Groups.size() == 0) {
				nameRep = findRep(namesInBin);
				newNamesOutput << i << '\t' << nameRep << '\t' << processList->get(i) << endl;
			}else{
				map<string, vector<string> > NamesInGroup;
				for (int j=0; j<Groups.size(); j++) { //initialize groups
					NamesInGroup[Groups[j]].resize(0);
				}
				
				for (int j=0; j<namesInBin.size(); j++) {
					string thisgroup = groupMap->getGroup(namesInBin[j]);
					
					if (thisgroup == "not found") { m->mothurOut(namesInBin[j] + " is not in your groupfile, please correct."); m->mothurOutEndLine(); m->control_pressed = true; }
					
					if (m->inUsersGroups(thisgroup, Groups)) { //add this name to correct group
						NamesInGroup[thisgroup].push_back(namesInBin[j]);
					}
				}
				
				//get rep for each group in otu
				for (int j=0; j<Groups.size(); j++) {
					if (NamesInGroup[Groups[j]].size() != 0) { //are there members from this group in this otu?
						//get rep for each group
						nameRep = findRep(NamesInGroup[Groups[j]]);
						
						//output group rep and other members of this group
						(*(filehandles[Groups[j]])) << i << '\t' << nameRep << '\t';
						
						for (int k=0; k<NamesInGroup[Groups[j]].size()-1; k++) {//output list of names in this otu from this group
							(*(filehandles[Groups[j]])) << NamesInGroup[Groups[j]][k] << ",";
						}
						//output last name
						(*(filehandles[Groups[j]])) << NamesInGroup[Groups[j]][NamesInGroup[Groups[j]].size()-1] << endl;
					}
				}
			}
		}
		
		if (Groups.size() == 0) { //you don't want to use groups
			newNamesOutput.close();
		}else{
			for (int i=0; i<Groups.size(); i++) {
				(*(filehandles[Groups[i]])).close();
				delete filehandles[Groups[i]];
			}
		}
		
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "GetOTURepCommand", "process");
		exit(1);
	}
}
//**********************************************************************************************************************
int GetOTURepCommand::processNames(string filename, string label) {
	try{

		//create output file
		if (outputDir == "") { outputDir += m->hasPath(listfile); }
		string outputFileName = outputDir + m->getRootName(m->getSimpleName(listfile)) + label + ".rep.fasta";
		m->openOutputFile(outputFileName, out);
		vector<repStruct> reps;
		outputNames.push_back(outputFileName); outputTypes["fasta"].push_back(outputFileName);
		
		ofstream out2;
		string tempNameFile = filename + ".temp";
		m->openOutputFile(tempNameFile, out2);
		
		ifstream in;
		m->openInputFile(filename, in);
		
		int i = 0;
		while (!in.eof()) {
			string rep, binnames;
			in >> i >> rep >> binnames; m->gobble(in);
			out2 << rep << '\t' << binnames << endl;
			
			vector<string> names;
			m->splitAtComma(binnames, names);
			int binsize = names.size();
			
			//if you have a groupfile
			string group = "";
			if (groupfile != "") {
				map<string, string> groups;
				map<string, string>::iterator groupIt;
				
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

			
			//print out name and sequence for that bin
			string sequence = fasta->getSequence(rep);

			if (sequence != "not found") {
				if (sorted == "") { //print them out
					rep = rep + "\t" + toString(i+1);
					rep = rep + "|" + toString(binsize);
					if (groupfile != "") {
						rep = rep + "|" + group;
					}
					out << ">" << rep << endl;
					out << sequence << endl;
				}else { //save them
					repStruct newRep(rep, i+1, binsize, group);
					reps.push_back(newRep);
				}
			}else { 
				m->mothurOut(rep + " is missing from your fasta or name file, ignoring. Please correct."); m->mothurOutEndLine(); 
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
				string outputName = reps[i].name + "\t" + toString(reps[i].bin);
				outputName = outputName + "|" + toString(reps[i].size);
				if (groupfile != "") {
					outputName = outputName + "|" + reps[i].group;
				}
				out << ">" << outputName << endl;
				out << sequence << endl;
			}
		}
		
		in.close();
		out.close();
		out2.close();
		
		m->mothurRemove(filename);
		rename(tempNameFile.c_str(), filename.c_str());
		
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "GetOTURepCommand", "processNames");
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

