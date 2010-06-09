/*
 *  splitabundcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 5/17/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "splitabundcommand.h"

//**********************************************************************************************************************
SplitAbundCommand::SplitAbundCommand(string option)  {
	try {
		abort = false;
		allLines = 1;
			
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"name","group","label","accnos","groups","fasta","cutoff","outputdir","inputdir"}; //"list",
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
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
				
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}

			}

			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";	}

			//check for required parameters
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { abort = true; }
			else if (listfile == "not found") { listfile = ""; }
			else{ inputFile = listfile; }	
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }
			else if (namefile == "not found") { namefile = ""; }	
			else{ inputFile = namefile; }	
		
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") { fastafile = ""; m->mothurOut("fasta is a required parameter for the split.abund command. "); m->mothurOutEndLine(); abort = true;  }	
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") {  groupfile = ""; abort = true; }	
			else if (groupfile == "not found") { groupfile = ""; }
			else {  
				groupMap = new GroupMap(groupfile);
				
				int error = groupMap->readMap();
				if (error == 1) { abort = true; }
	
			}
			
			groups = validParameter.validFile(parameters, "groups", false);		
			if (groups == "not found") { groups = ""; }
			else if (groups == "all") { 
				if (groupfile != "") {  Groups = groupMap->namesOfGroups;  } 
				else {  m->mothurOut("You cannot select groups without a valid groupfile, I will disregard your groups selection. "); m->mothurOutEndLine(); groups = "";   }
			}else { 
				splitAtDash(groups, Groups);
			}
			
			if ((groupfile == "") && (groups != "")) {  m->mothurOut("You cannot select groups without a valid groupfile, I will disregard your groups selection. "); m->mothurOutEndLine(); groups = "";  Groups.clear(); }
			
			//do you have all files needed
			if ((listfile == "") && (namefile == "")) { m->mothurOut("You must either a listfile or a namefile for the split.abund command. "); m->mothurOutEndLine(); abort = true;  }
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.validFile(parameters, "label", false);			
			if (label == "not found") { label = "";  allLines = 1; }
			else { 
				if(label != "all") {  splitAtDash(label, labels);  allLines = 0;  }
				else { allLines = 1;  }
			}
			
			string temp = validParameter.validFile(parameters, "accnos", false);		if (temp == "not found") { temp = "F"; }
			accnos = isTrue(temp); 
			
			temp = validParameter.validFile(parameters, "cutoff", false);				if (temp == "not found") { temp = "0"; }
			convert(temp, cutoff); 

			if (cutoff == 0) {  m->mothurOut("You must provide a cutoff to qualify what is abundant for the split.abund command. "); m->mothurOutEndLine(); abort = true;  }

		}

	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "SplitAbundCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
void SplitAbundCommand::help(){
	try {
		m->mothurOut("The split.abund command reads a fasta file and a list or a names file splits the sequences into rare and abundant groups. \n");
		m->mothurOut("The split.abund command parameters are list, name, cutoff, group, label, groups and accnos.\n");
		m->mothurOut("The fasta and a list or name parameter are required, and you must provide a cutoff value.\n");
		m->mothurOut("The cutoff parameter is used to qualify what is abundant and rare.\n");
		m->mothurOut("The group parameter allows you to parse a group file into rare and abundant groups.\n");
		m->mothurOut("The label parameter is used to read specific labels in your listfile you want to use.\n");
		m->mothurOut("The accnos parameter allows you to output a .rare.accnos and .abund.accnos files to use with the get.seqs and remove.seqs commands.\n");
		m->mothurOut("The groups parameter allows you to parse the files into rare and abundant files by group.  \n");
		m->mothurOut("For example if you set groups=A-B-C, you will get a .A.abund, .A.rare, .B.abund, .B.rare, .C.abund, .C.rare files.  \n");
		m->mothurOut("If you want .abund and .rare files for all groups, set groups=all.  \n");
		m->mothurOut("The split.abund command should be used in the following format: split.abund(list=yourListFile, group=yourGroupFile, label=yourLabels, cutoff=yourCutoff).\n");
		m->mothurOut("Example: split.abundt(list=abrecovery.fn.list, group=abrecovery.groups, label=0.03, cutoff=2).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. list), '=' and parameters (i.e.yourListfile).\n\n");

	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "help");
		exit(1);
	}
}
//**********************************************************************************************************************
SplitAbundCommand::~SplitAbundCommand(){ 
	if (groupfile != "") {  delete groupMap;  } 
}
//**********************************************************************************************************************
int SplitAbundCommand::execute(){
	try {
	
		if (abort == true) {	return 0;	}
		
		if (listfile != "") { //you are using a listfile to determine abundance
		
			//remove old files so you can append later....
			string fileroot = outputDir + getRootName(getSimpleName(listfile));
			if (Groups.size() == 0) {
				remove((fileroot + "rare.list").c_str());
				remove((fileroot + "abund.list").c_str());
				
				wroteListFile["rare"] = false;
				wroteListFile["abund"] = false;
			}else{
				for (int i=0; i<Groups.size(); i++) {
					remove((fileroot + Groups[i] + ".rare.list").c_str());
					remove((fileroot + Groups[i] + ".abund.list").c_str());
			
					wroteListFile[(Groups[i] + ".rare")] = false;
					wroteListFile[(Groups[i] + ".abund")] = false;
				}
			}
			
			//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
			set<string> processedLabels;
			set<string> userLabels = labels;	
			
			input = new InputData(listfile, "list");
			list = input->getListVector();
			string lastLabel = list->getLabel();
			
			//do you have a namefile or do we need to similate one?
			if (namefile != "") {  readNamesFile();		}
			else				{ createNameMap(list);	}
			
			if (m->control_pressed) { delete input; delete list; for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }
			
			while((list != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			
				if (m->control_pressed) { delete input; delete list; for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }
				
				if(allLines == 1 || labels.count(list->getLabel()) == 1){
						
						m->mothurOut(list->getLabel()); m->mothurOutEndLine();
						splitList(list);
											
						processedLabels.insert(list->getLabel());
						userLabels.erase(list->getLabel());
				}
				
				if ((anyLabelsToProcess(list->getLabel(), userLabels, "") == true) && (processedLabels.count(lastLabel) != 1)) {
						string saveLabel = list->getLabel();
						
						delete list;
						list = input->getListVector(lastLabel); //get new list vector to process
						
						m->mothurOut(list->getLabel()); m->mothurOutEndLine();
						splitList(list);
						
						processedLabels.insert(list->getLabel());
						userLabels.erase(list->getLabel());
						
						//restore real lastlabel to save below
						list->setLabel(saveLabel);
				}
				
			
				lastLabel = list->getLabel();
					
				delete list;
				list = input->getListVector(); //get new list vector to process
			}
			
			if (m->control_pressed) { delete input;  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }
			
			//output error messages about any remaining user labels
			set<string>::iterator it;
			bool needToRun = false;
			for (it = userLabels.begin(); it != userLabels.end(); it++) {  
				m->mothurOut("Your file does not include the label " + *it); 
				if (processedLabels.count(lastLabel) != 1) {
					m->mothurOut(". I will use " + lastLabel + "."); m->mothurOutEndLine();
					needToRun = true;
				}else {
					m->mothurOut(". Please refer to " + lastLabel + "."); m->mothurOutEndLine();
				}

			}
			
			if (m->control_pressed) { delete input;  for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }
			
			//run last label if you need to
			if (needToRun == true)  {
				if (list != NULL) {	delete list;	}
				list = input->getListVector(lastLabel); //get new list vector to process
				
				m->mothurOut(list->getLabel()); m->mothurOutEndLine();
				splitList(list);		
				
				delete list;
			}
			
			delete input;
			
			for (map<string, bool>::iterator itBool = wroteListFile.begin(); itBool != wroteListFile.end(); itBool++) {
				string filename = fileroot + itBool->first;
				if (itBool->second) { //we wrote to this file
					outputNames.push_back(filename);
				}else{
					remove(filename.c_str());
				}
			}
			
			if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); }	return 0;	}

									
		}else { //you are using the namefile to determine abundance

			splitNames(); 
			writeNames();
			
			string tag = "";
			if (groupfile != "")				{  parseGroup(tag);		}
			if (accnos)							{  writeAccnos(tag);	}
			if (fastafile != "")				{  parseFasta(tag);		}
		}

		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "execute");
		exit(1);
	}
}
/**********************************************************************************************************************/
int SplitAbundCommand::splitList(ListVector* thisList) {
	try {
		rareNames.clear();
		abundNames.clear();
		
		//get rareNames and abundNames
		for (int i = 0; i < thisList->getNumBins(); i++) {
			if (m->control_pressed) { return 0; }
			
			string bin = thisList->get(i);
						
			vector<string> names;
			splitAtComma(bin, names);  //parses bin into individual sequence names
			int size = names.size();
				
			if (size <= cutoff) {
				for (int j = 0; j < names.size(); j++) {  rareNames.insert(names[j]);  }
			}else{
				for (int j = 0; j < names.size(); j++) {  abundNames.insert(names[j]);  }
			}
		}//end for

		writeList(thisList);
		
		string tag = thisList->getLabel() + ".";
		if (groupfile != "")				{  parseGroup(tag);		}
		if (accnos)							{  writeAccnos(tag);	}
		if (fastafile != "")				{  parseFasta(tag);		}
		
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "splitList");
		exit(1);
	}
}
/**********************************************************************************************************************/
int SplitAbundCommand::writeList(ListVector* thisList) { 
	try {
		
		map<string, ofstream*> filehandles;
		
		if (Groups.size() == 0) {
			SAbundVector* sabund = new SAbundVector();
			*sabund = thisList->getSAbundVector();
		
			//find out how many bins are rare and how many are abundant so you can process the list vector one bin at a time
			// and don't have to store the bins until you are done with the whole vector, this save alot of space.
			int numRareBins = 0;
			for (int i = 0; i <= sabund->getMaxRank(); i++) {
				if (i > cutoff) { break; }
				numRareBins += sabund->get(i);
			}
			int numAbundBins = thisList->getNumBins() - numRareBins;
			delete sabund;

			ofstream aout;
			ofstream rout;
			
			if (rareNames.size() != 0) {
				string rare = outputDir + getRootName(getSimpleName(listfile))  + ".rare.list";
				wroteListFile["rare"] = true;
				openOutputFileAppend(rare, rout);
				rout << thisList->getLabel() << '\t' << numRareBins << '\t';
			}
			
			if (abundNames.size() != 0) {
				string abund = outputDir + getRootName(getSimpleName(listfile))  + ".abund.list";
				wroteListFile["abund"] = true;
				openOutputFileAppend(abund, aout);
				rout << thisList->getLabel() << '\t' << numAbundBins << '\t';
			}

			for (int i = 0; i < thisList->getNumBins(); i++) {
				if (m->control_pressed) { break; }
			
				string bin = list->get(i); 
			
				int size = getNumNames(bin);
			
				if (size <= cutoff) {  rout << bin << '\t';  }
				else				{  aout << bin << '\t'; }
			}
			
			if (rareNames.size() != 0)	{ rout << endl; rout.close(); }
			if (abundNames.size() != 0) { aout << endl; aout.close(); }

		}else{ //parse names by abundance and group
			string fileroot =  outputDir + getRootName(getSimpleName(listfile));
			ofstream* temp;
			ofstream* temp2;
			map<string, bool> wroteFile;
			map<string, ofstream*> filehandles;
			map<string, ofstream*>::iterator it3;

			for (int i=0; i<Groups.size(); i++) {
				temp = new ofstream;
				filehandles[Groups[i]+".rare"] = temp;
				temp2 = new ofstream;
				filehandles[Groups[i]+".abund"] = temp2;
				
				openOutputFileAppend(fileroot + Groups[i] + ".rare.list", *(filehandles[Groups[i]+".rare"]));
				openOutputFileAppend(fileroot + Groups[i] + ".abund.list", *(filehandles[Groups[i]+".abund"]));
			}
			
			map<string, string> groupVector;
			map<string, string>::iterator itGroup;
			map<string, int> groupNumBins;
		
			for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {
				groupNumBins[it3->first] = 0;
				groupVector[it3->first] = "";
			}
		
			for (int i = 0; i < thisList->getNumBins(); i++) {
				if (m->control_pressed) { break; }
			
				map<string, string> groupBins;
				string bin = list->get(i); 
			
				vector<string> names;
				splitAtComma(bin, names);  //parses bin into individual sequence names
			
				//parse bin into list of sequences in each group
				for (int j = 0; j < names.size(); j++) {
					string rareAbund;
					if (rareNames.count(names[j]) != 0) { //you are a rare name
						rareAbund = ".rare";
					}else{ //you are a abund name
						rareAbund = ".abund";
					}
					
					string group = groupMap->getGroup(names[j]);
				
					if (inUsersGroups(group, Groups)) { //only add if this is in a group we want
						itGroup = groupBins.find(group+rareAbund);
						if(itGroup == groupBins.end()) {
							groupBins[group+rareAbund] = names[j];  //add first name
							groupNumBins[group+rareAbund]++;
						}else{ //add another name
							groupBins[group+rareAbund] +=  "," + names[j];
						}
					}else if(group == "not found") {
						m->mothurOut(names[j] + " is not in your groupfile. Ignoring."); m->mothurOutEndLine();
					}
				}
			
			
				for (itGroup = groupBins.begin(); itGroup != groupBins.end(); itGroup++) {
					groupVector[itGroup->first] +=  itGroup->second + '\t'; 
				}
			}
			
			//end list vector
			for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) {
				(*(filehandles[it3->first])) << thisList->getLabel() << '\t' << groupNumBins[it3->first] << '\t' << groupVector[it3->first] << endl;  // label numBins  listvector for that group
				wroteListFile[it3->first] = true;
				(*(filehandles[it3->first])).close();
				delete it3->second;
			}
		}
		
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "writeList");
		exit(1);
	}
}
/**********************************************************************************************************************/
int SplitAbundCommand::splitNames() { //namefile
	try {
		
		rareNames.clear();
		abundNames.clear();	
			
		//open input file
		ifstream in;
		openInputFile(namefile, in);
		
		while (!in.eof()) {
			if (m->control_pressed) { break; }
			
			string firstCol, secondCol;
			in >> firstCol >> secondCol; gobble(in);
			
			nameMap[firstCol] = secondCol;
			
			int size = getNumNames(secondCol);
				
			if (size <= cutoff) {
				rareNames.insert(firstCol); 
			}else{
				abundNames.insert(firstCol); 
			}
		}
		in.close();
				
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "splitNames");
		exit(1);
	}
}
/**********************************************************************************************************************/
int SplitAbundCommand::readNamesFile() { 
	try {
		//open input file
		ifstream in;
		openInputFile(namefile, in);
		
		while (!in.eof()) {
			if (m->control_pressed) { break; }
			
			string firstCol, secondCol;
			in >> firstCol >> secondCol; gobble(in);
			
			nameMap[firstCol] = secondCol;
		}
		in.close();
				
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "readNamesFile");
		exit(1);
	}
}
/**********************************************************************************************************************/
int SplitAbundCommand::createNameMap(ListVector* thisList) {
	try {
		
		if (thisList != NULL) {
			for (int i = 0; i < thisList->getNumBins(); i++) {
				if (m->control_pressed) { return 0; }
				
				string bin = thisList->get(i);
							
				vector<string> names;
				splitAtComma(bin, names);  //parses bin into individual sequence names
				
				for (int j = 0; j < names.size(); j++) {  nameMap[names[j]] = names[j];  }
			}//end for
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "createNameMap");
		exit(1);
	}
}
/**********************************************************************************************************************/
int SplitAbundCommand::writeNames() { //namefile
	try {
		
		map<string, ofstream*> filehandles;

		if (Groups.size() == 0) {
			ofstream aout;
			ofstream rout;
			
			if (rareNames.size() != 0) {
				string rare = outputDir + getRootName(getSimpleName(namefile))  + "rare.names";
				openOutputFile(rare, rout);
				outputNames.push_back(rare);
				
				for (set<string>::iterator itRare = rareNames.begin(); itRare != rareNames.end(); itRare++) {
					rout << (*itRare) << '\t' << nameMap[(*itRare)] << endl;
				}
				rout.close();
			}
			
			if (abundNames.size() != 0) {
				string abund = outputDir + getRootName(getSimpleName(namefile))  + "abund.names";
				openOutputFile(abund, aout);
				outputNames.push_back(abund);
				
				for (set<string>::iterator itAbund = abundNames.begin(); itAbund != abundNames.end(); itAbund++) {
					aout << (*itAbund) << '\t' << nameMap[(*itAbund)] << endl;
				}
				aout.close();
			}

		}else{ //parse names by abundance and group
			string fileroot =  outputDir + getRootName(getSimpleName(namefile));
			ofstream* temp;
			ofstream* temp2;
			map<string, bool> wroteFile;
			map<string, ofstream*> filehandles;
			map<string, ofstream*>::iterator it3;

			for (int i=0; i<Groups.size(); i++) {
				temp = new ofstream;
				filehandles[Groups[i]+".rare"] = temp;
				temp2 = new ofstream;
				filehandles[Groups[i]+".abund"] = temp2;
				
				openOutputFile(fileroot + Groups[i] + ".rare.names", *(filehandles[Groups[i]+".rare"]));
				openOutputFile(fileroot + Groups[i] + ".abund.names", *(filehandles[Groups[i]+".abund"]));
				
				wroteFile[Groups[i] + ".rare"] = false;
				wroteFile[Groups[i] + ".abund"] = false;
			}
			
			for (map<string, string>::iterator itName = nameMap.begin(); itName != nameMap.end(); itName++) {				
				vector<string> names;
				splitAtComma(itName->second, names);  //parses bin into individual sequence names
				
				string rareAbund;
				if (rareNames.count(itName->first) != 0) { //you are a rare name
						rareAbund = ".rare";
				}else{ //you are a abund name
						rareAbund = ".abund";
				}
				
				map<string, string> outputStrings;
				map<string, string>::iterator itout;
				for (int i = 0; i < names.size(); i++) {
					
					string group = groupMap->getGroup(names[i]);
					
					if (inUsersGroups(group, Groups)) { //only add if this is in a group we want
						itout = outputStrings.find(group+rareAbund);
						if (itout == outputStrings.end()) {  
							outputStrings[group+rareAbund] = names[i] + '\t' + names[i];
						}else {   outputStrings[group+rareAbund] += "," + names[i]; }
					}else if(group == "not found") {
						m->mothurOut(names[i] + " is not in your groupfile. Ignoring."); m->mothurOutEndLine();
					}
				}
				
				for (itout = outputStrings.begin(); itout != outputStrings.end(); itout++) { 
					*(filehandles[itout->first]) << itout->second << endl;
					wroteFile[itout->first] = true;
				}
			}
			
			
			for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) { 
				(*(filehandles[it3->first])).close();
				if (wroteFile[it3->first] == true) {  outputNames.push_back(fileroot + it3->first + ".names");  }
				else { remove((it3->first).c_str()); }
				delete it3->second;
			}
		}
				
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "writeNames");
		exit(1);
	}
}
/**********************************************************************************************************************/
//just write the unique names - if a namesfile is given
int SplitAbundCommand::writeAccnos(string tag) { 
	try {
		
		map<string, ofstream*> filehandles;
		
		if (Groups.size() == 0) {
			ofstream aout;
			ofstream rout;
			
			if (rareNames.size() != 0) {
				string rare = outputDir + getRootName(getSimpleName(inputFile))  + tag + "rare.accnos";
				openOutputFile(rare, rout);
				outputNames.push_back(rare);
				
				for (set<string>::iterator itRare = rareNames.begin(); itRare != rareNames.end(); itRare++) {
					rout << (*itRare) << endl;
				}
				rout.close();
			}
			
			if (abundNames.size() != 0) {
				string abund = outputDir + getRootName(getSimpleName(inputFile)) + tag  + "abund.accnos";
				openOutputFile(abund, aout);
				outputNames.push_back(abund);
				
				for (set<string>::iterator itAbund = abundNames.begin(); itAbund != abundNames.end(); itAbund++) {
					aout << (*itAbund) << endl;
				}
				aout.close();
			}
		}else{ //parse names by abundance and group
			string fileroot =  outputDir + getRootName(getSimpleName(inputFile));
			ofstream* temp;
			ofstream* temp2;
			map<string, bool> wroteFile;
			map<string, ofstream*> filehandles;
			map<string, ofstream*>::iterator it3;
			
			for (int i=0; i<Groups.size(); i++) {
				temp = new ofstream;
				filehandles[Groups[i]+".rare"] = temp;
				temp2 = new ofstream;
				filehandles[Groups[i]+".abund"] = temp2;
				
				openOutputFile(fileroot + tag + Groups[i] + ".rare.accnos", *(filehandles[Groups[i]+".rare"]));
				openOutputFile(fileroot + tag + Groups[i] + ".abund.accnos", *(filehandles[Groups[i]+".abund"]));
				
				wroteFile[Groups[i] + ".rare"] = false;
				wroteFile[Groups[i] + ".abund"] = false;
			}
			
			//write rare
			for (set<string>::iterator itRare = rareNames.begin(); itRare != rareNames.end(); itRare++) {
					string group = groupMap->getGroup(*itRare);
					
					if (inUsersGroups(group, Groups)) { //only add if this is in a group we want
						*(filehandles[group+".rare"]) << *itRare << endl;
						wroteFile[group+".rare"] = true;
					}
			}
				
			//write abund	
			for (set<string>::iterator itAbund = abundNames.begin(); itAbund != abundNames.end(); itAbund++) {
					string group = groupMap->getGroup(*itAbund);
					
					if (inUsersGroups(group, Groups)) { //only add if this is in a group we want
						*(filehandles[group+".abund"]) << *itAbund << endl;
						wroteFile[group+".abund"] = true;
					}
			}
			
			//close files
			for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) { 
				(*(filehandles[it3->first])).close();
				if (wroteFile[it3->first] == true) {  outputNames.push_back(fileroot + tag + it3->first + ".accnos");  }
				else { remove((fileroot + tag + it3->first + ".accnos").c_str()); }
				delete it3->second;
			}
		}
				
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "writeAccnos");
		exit(1);
	}
}
/**********************************************************************************************************************/
int SplitAbundCommand::parseGroup(string tag) { //namefile
	try {
		
		map<string, ofstream*> filehandles;
	
		if (Groups.size() == 0) {
			ofstream aout;
			ofstream rout;
			
			if (rareNames.size() != 0) {
				string rare = outputDir + getRootName(getSimpleName(groupfile))  + tag + "rare.groups";
				openOutputFile(rare, rout);
				outputNames.push_back(rare);
			}
			
			if (abundNames.size() != 0) {
				string abund = outputDir + getRootName(getSimpleName(groupfile))  + tag + "abund.groups";
				openOutputFile(abund, aout);
				outputNames.push_back(abund);
			}

				
			for (map<string, string>::iterator itName = nameMap.begin(); itName != nameMap.end(); itName++) {				
				vector<string> names;
				splitAtComma(itName->second, names);  //parses bin into individual sequence names
				
				for (int i = 0; i < names.size(); i++) {
				
					string group = groupMap->getGroup(names[i]);
				
					if (group == "not found") { 
						m->mothurOut(names[i] + " is not in your groupfile, ignoring, please correct."); m->mothurOutEndLine();
					}else {
						if (rareNames.count(itName->first) != 0) { //you are a rare name
							rout << names[i] << '\t' << group << endl;
						}else{ //you are a abund name
							rout << names[i] << '\t' << group << endl;
						}
					}
				}
			}
			
			if (rareNames.size() != 0)	{ rout.close(); }
			if (abundNames.size() != 0) { aout.close(); }

		}else{ //parse names by abundance and group
			string fileroot =  outputDir + getRootName(getSimpleName(groupfile));
			ofstream* temp;
			ofstream* temp2;
			map<string, bool> wroteFile;
			map<string, ofstream*> filehandles;
			map<string, ofstream*>::iterator it3;

			for (int i=0; i<Groups.size(); i++) {
				temp = new ofstream;
				filehandles[Groups[i]+".rare"] = temp;
				temp2 = new ofstream;
				filehandles[Groups[i]+".abund"] = temp2;
				
				openOutputFile(fileroot + tag + Groups[i] + ".rare.groups", *(filehandles[Groups[i]+".rare"]));
				openOutputFile(fileroot + tag + Groups[i] + ".abund.groups", *(filehandles[Groups[i]+".abund"]));
				
				wroteFile[Groups[i] + ".rare"] = false;
				wroteFile[Groups[i] + ".abund"] = false;
			}
			
			for (map<string, string>::iterator itName = nameMap.begin(); itName != nameMap.end(); itName++) {				
				vector<string> names;
				splitAtComma(itName->second, names);  //parses bin into individual sequence names
				
				string rareAbund;
				if (rareNames.count(itName->first) != 0) { //you are a rare name
					rareAbund = ".rare";
				}else{ //you are a abund name
					rareAbund = ".abund";
				}
				
				for (int i = 0; i < names.size(); i++) {
				
					string group = groupMap->getGroup(names[i]);
									
					if (inUsersGroups(group, Groups)) { //only add if this is in a group we want
						*(filehandles[group+rareAbund]) << names[i] << '\t' << group << endl;
						wroteFile[group+rareAbund] = true;
					}
				}
			}
			
			for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) { 
				(*(filehandles[it3->first])).close();
				if (wroteFile[it3->first] == true) {  outputNames.push_back(fileroot + tag + it3->first + ".groups");  }
				else { remove((fileroot + tag + it3->first + ".groups").c_str()); }
				delete it3->second;
			}
		}
				
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "parseGroups");
		exit(1);
	}
}
/**********************************************************************************************************************/
int SplitAbundCommand::parseFasta(string tag) { //namefile
	try {
		
		map<string, ofstream*> filehandles;
		
		if (Groups.size() == 0) {
			ofstream aout;
			ofstream rout;
			
			if (rareNames.size() != 0) {
				string rare = outputDir + getRootName(getSimpleName(fastafile))  + tag + "rare.fasta";
				openOutputFile(rare, rout);
				outputNames.push_back(rare);
			}
			
			if (abundNames.size() != 0) {
				string abund = outputDir + getRootName(getSimpleName(fastafile))  + tag + "abund.fasta";
				openOutputFile(abund, aout);
				outputNames.push_back(abund);
			}

				
			//open input file
			ifstream in;
			openInputFile(fastafile, in);
	
			while (!in.eof()) {
				if (m->control_pressed) { break; }
		
				Sequence seq(in); gobble(in);
				
				if (seq.getName() != "") { 
					
					map<string, string>::iterator itNames;
					
					itNames = nameMap.find(seq.getName());
					
					if (itNames == nameMap.end()) {
						m->mothurOut(seq.getName() + " is not in your namesfile, ignoring."); m->mothurOutEndLine();
					}else{
						if (rareNames.count(seq.getName()) != 0) { //you are a rare name
							seq.printSequence(rout);
						}else{ //you are a abund name
							seq.printSequence(aout);
						}
					}
				}
			}
			in.close();
			if (rareNames.size() != 0)	{ rout.close(); }
			if (abundNames.size() != 0) { aout.close(); }

		}else{ //parse names by abundance and group
			string fileroot =  outputDir + getRootName(getSimpleName(fastafile));
			ofstream* temp;
			ofstream* temp2;
			map<string, bool> wroteFile;
			map<string, ofstream*> filehandles;
			map<string, ofstream*>::iterator it3;

			for (int i=0; i<Groups.size(); i++) {
				temp = new ofstream;
				filehandles[Groups[i]+".rare"] = temp;
				temp2 = new ofstream;
				filehandles[Groups[i]+".abund"] = temp2;
				
				openOutputFile(fileroot + tag + Groups[i] + ".rare.fasta", *(filehandles[Groups[i]+".rare"]));
				openOutputFile(fileroot + tag + Groups[i] + ".abund.fasta", *(filehandles[Groups[i]+".abund"]));
				
				wroteFile[Groups[i] + ".rare"] = false;
				wroteFile[Groups[i] + ".abund"] = false;
			}
			
			//open input file
			ifstream in;
			openInputFile(fastafile, in);
	
			while (!in.eof()) {
				if (m->control_pressed) { break; }
		
				Sequence seq(in); gobble(in);
				
				if (seq.getName() != "") { 
					map<string, string>::iterator itNames = nameMap.find(seq.getName());
					
					if (itNames == nameMap.end()) {
						m->mothurOut(seq.getName() + " is not in your namesfile, ignoring."); m->mothurOutEndLine();
					}else{
						vector<string> names;
						splitAtComma(itNames->second, names);  //parses bin into individual sequence names
				
						string rareAbund;
						if (rareNames.count(itNames->first) != 0) { //you are a rare name
							rareAbund = ".rare";
						}else{ //you are a abund name
							rareAbund = ".abund";
						}
				
						for (int i = 0; i < names.size(); i++) {
				
							string group = groupMap->getGroup(seq.getName());
					
							if (inUsersGroups(group, Groups)) { //only add if this is in a group we want
								seq.printSequence(*(filehandles[group+rareAbund]));
								wroteFile[group+rareAbund] = true;
							}else if(group == "not found") {
								m->mothurOut(seq.getName() + " is not in your groupfile. Ignoring."); m->mothurOutEndLine();
							}
						}
					}
				}
			}
			in.close();
			
			for (it3 = filehandles.begin(); it3 != filehandles.end(); it3++) { 
				(*(filehandles[it3->first])).close();
				if (wroteFile[it3->first] == true) {  outputNames.push_back(fileroot + tag + it3->first + ".fasta");  }
				else { remove((fileroot + tag + it3->first + ".fasta").c_str()); }
				delete it3->second;
			}
		}
				
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "parseFasta");
		exit(1);
	}
}
/**********************************************************************************************************************/

