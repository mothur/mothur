/*
 *  removeseqscommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "removeseqscommand.h"
#include "sequence.hpp"
#include "listvector.hpp"

//**********************************************************************************************************************
vector<string> RemoveSeqsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(pname);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(pgroup);
		CommandParameter plist("list", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(plist);
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(ptaxonomy);
		CommandParameter palignreport("alignreport", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(palignreport);
		CommandParameter pqfile("qfile", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(pqfile);
		CommandParameter paccnos("accnos", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(paccnos);
		CommandParameter pdups("dups", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(pdups);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string RemoveSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The remove.seqs command reads an .accnos file and at least one of the following file types: fasta, name, group, list, taxonomy, quality or alignreport file.\n";
		helpString += "It outputs a file containing the sequences NOT in the .accnos file.\n";
		helpString += "The remove.seqs command parameters are accnos, fasta, name, group, list, taxonomy, qfile, alignreport and dups.  You must provide accnos and at least one of the file parameters.\n";
		helpString += "The dups parameter allows you to remove the entire line from a name file if you remove any name from the line. default=true. \n";
		helpString += "The remove.seqs command should be in the following format: remove.seqs(accnos=yourAccnos, fasta=yourFasta).\n";
		helpString += "Example remove.seqs(accnos=amazon.accnos, fasta=amazon.fasta).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "getHelpString");
		exit(1);
	}
}


//**********************************************************************************************************************
RemoveSeqsCommand::RemoveSeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["taxonomy"] = tempOutNames;
		outputTypes["name"] = tempOutNames;
		outputTypes["group"] = tempOutNames;
		outputTypes["alignreport"] = tempOutNames;
		outputTypes["list"] = tempOutNames;
		outputTypes["qfile"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "RemoveSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
RemoveSeqsCommand::RemoveSeqsCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
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
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["taxonomy"] = tempOutNames;
			outputTypes["name"] = tempOutNames;
			outputTypes["group"] = tempOutNames;
			outputTypes["alignreport"] = tempOutNames;
			outputTypes["list"] = tempOutNames;
			outputTypes["qfile"] = tempOutNames;
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("alignreport");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["alignreport"] = inputDir + it->second;		}
				}
				
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("accnos");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["accnos"] = inputDir + it->second;		}
				}
				
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
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
				
				it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
				}
				
				it = parameters.find("qfile");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["qfile"] = inputDir + it->second;		}
				}
			}

			
			//check for required parameters
			accnosfile = validParameter.validFile(parameters, "accnos", true);
			if (accnosfile == "not open") { abort = true; }
			else if (accnosfile == "not found") {  
				accnosfile = m->getAccnosFile(); 
				if (accnosfile != "") {  m->mothurOut("Using " + accnosfile + " as input file for the accnos parameter."); m->mothurOutEndLine(); }
				else { 
					m->mothurOut("You have no valid accnos file and accnos is required."); m->mothurOutEndLine(); 
					abort = true;
				}  
			}	
			
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") {  fastafile = "";  }	
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { abort = true; }
			else if (namefile == "not found") {  namefile = "";  }	
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }
			else if (groupfile == "not found") {  groupfile = "";  }	
			
			alignfile = validParameter.validFile(parameters, "alignreport", true);
			if (alignfile == "not open") { abort = true; }
			else if (alignfile == "not found") {  alignfile = "";  }
			
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { abort = true; }
			else if (listfile == "not found") {  listfile = "";  }
			
			taxfile = validParameter.validFile(parameters, "taxonomy", true);
			if (taxfile == "not open") { abort = true; }
			else if (taxfile == "not found") {  taxfile = "";  }
			
			qualfile = validParameter.validFile(parameters, "qfile", true);
			if (qualfile == "not open") { abort = true; }
			else if (qualfile == "not found") {  qualfile = "";  }			

			
			string usedDups = "true";
			string temp = validParameter.validFile(parameters, "dups", false);	
			if (temp == "not found") { 
				if (namefile != "") {  temp = "true";					}
				else				{  temp = "false"; usedDups = "";	}
			}
			dups = m->isTrue(temp);
			
			if ((fastafile == "") && (namefile == "") && (groupfile == "") && (alignfile == "") && (listfile == "") && (taxfile == "") && (qualfile == ""))  { m->mothurOut("You must provide at least one of the following: fasta, name, group, taxonomy, quality, alignreport or list."); m->mothurOutEndLine(); abort = true; }
			
		}

	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "RemoveSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int RemoveSeqsCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//get names you want to keep
		readAccnos();
		
		if (m->control_pressed) { return 0; }
		
		//read through the correct file and output lines you want to keep
		if (namefile != "")			{		readName();		}
		if (fastafile != "")		{		readFasta();	}
		if (groupfile != "")		{		readGroup();	}
		if (alignfile != "")		{		readAlign();	}
		if (listfile != "")			{		readList();		}
		if (taxfile != "")			{		readTax();		}
		if (qualfile != "")			{		readQual();		}
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str()); } return 0; }
		
		m->mothurOut("Removed " + toString(names.size()) + " sequences."); m->mothurOutEndLine();
		
		if (outputNames.size() != 0) {
			m->mothurOutEndLine();
			m->mothurOut("Output File Names: "); m->mothurOutEndLine();
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
			m->mothurOutEndLine();
			
			//set fasta file as new current fastafile
			string current = "";
			itTypes = outputTypes.find("fasta");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
			}
			
			itTypes = outputTypes.find("name");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setNameFile(current); }
			}
			
			itTypes = outputTypes.find("group");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setGroupFile(current); }
			}
			
			itTypes = outputTypes.find("list");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setListFile(current); }
			}
			
			itTypes = outputTypes.find("taxonomy");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setTaxonomyFile(current); }
			}
			
			itTypes = outputTypes.find("qfile");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setQualFile(current); }
			}			
		}
		
		return 0;		
	}

	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
int RemoveSeqsCommand::readFasta(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(fastafile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(fastafile)) + "pick" + m->getExtension(fastafile);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(fastafile, in);
		string name;
		
		bool wroteSomething = false;
		
		while(!in.eof()){
			if (m->control_pressed) { in.close();  out.close();  remove(outputFileName.c_str());  return 0; }
			
			Sequence currSeq(in);
			name = currSeq.getName();
			
			if (name != "") {
				//if this name is in the accnos file
				if (names.count(name) == 0) {
					wroteSomething = true;
					
					currSeq.printSequence(out);
				}
			}
			m->gobble(in);
		}
		in.close();	
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the .accnos file."); m->mothurOutEndLine();  }
		outputTypes["fasta"].push_back(outputFileName);  outputNames.push_back(outputFileName);
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "readFasta");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveSeqsCommand::readQual(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(qualfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(qualfile)) + "pick" +  m->getExtension(qualfile);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		
		ifstream in;
		m->openInputFile(qualfile, in);
		string name;
		
		bool wroteSomething = false;
		
		
		while(!in.eof()){	
			string saveName = "";
			string name = "";
			string scores = "";
			
			in >> name; 
			
			if (name.length() != 0) { 
				saveName = name.substr(1);
				while (!in.eof())	{	
					char c = in.get(); 
					if (c == 10 || c == 13){	break;	}
					else { name += c; }	
				} 
				m->gobble(in);
			}
			
			while(in){
				char letter= in.get();
				if(letter == '>'){	in.putback(letter);	break;	}
				else{ scores += letter; }
			}
			
			m->gobble(in);
			
			if (names.count(saveName) == 0) {
				wroteSomething = true;
				
				out << name << endl << scores;
			}
			
			m->gobble(in);
		}
		in.close();
		out.close();
		
		
		if (wroteSomething == false) { m->mothurOut("Your file contains only sequences from the .accnos file."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName);  outputTypes["qfile"].push_back(outputFileName); 
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "readQual");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveSeqsCommand::readList(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(listfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(listfile)) + "pick" +  m->getExtension(listfile);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(listfile, in);
		
		bool wroteSomething = false;
		
		while(!in.eof()){
			//read in list vector
			ListVector list(in);
			
			//make a new list vector
			ListVector newList;
			newList.setLabel(list.getLabel());
			
			//for each bin
			for (int i = 0; i < list.getNumBins(); i++) {
				if (m->control_pressed) { in.close();  out.close();  remove(outputFileName.c_str());  return 0; }
			
				//parse out names that are in accnos file
				string binnames = list.get(i);
				
				string newNames = "";
				while (binnames.find_first_of(',') != -1) { 
					string name = binnames.substr(0,binnames.find_first_of(','));
					binnames = binnames.substr(binnames.find_first_of(',')+1, binnames.length());
					
					//if that name is in the .accnos file, add it
					if (names.count(name) == 0) {  newNames += name + ",";  }
				}
			
				//get last name
				if (names.count(binnames) == 0) {  newNames += binnames + ",";  }

				//if there are names in this bin add to new list
				if (newNames != "") {  
					newNames = newNames.substr(0, newNames.length()-1); //rip off extra comma
					newList.push_back(newNames);	
				}
			}
				
			//print new listvector
			if (newList.getNumBins() != 0) {
				wroteSomething = true;
				newList.print(out);
			}
			
			m->gobble(in);
		}
		in.close();	
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the .accnos file."); m->mothurOutEndLine();  }
		outputTypes["list"].push_back(outputFileName); outputNames.push_back(outputFileName);
				
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "readList");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveSeqsCommand::readName(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(namefile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(namefile)) + "pick" + m->getExtension(namefile);

		ofstream out;
		m->openOutputFile(outputFileName, out);

		ifstream in;
		m->openInputFile(namefile, in);
		string name, firstCol, secondCol;
		
		bool wroteSomething = false;
				
		
		while(!in.eof()){
			if (m->control_pressed) { in.close();  out.close();  remove(outputFileName.c_str());  return 0; }
			
			in >> firstCol;		m->gobble(in);		
			in >> secondCol;			
			
			vector<string> parsedNames;
			m->splitAtComma(secondCol, parsedNames);
			
			vector<string> validSecond;  validSecond.clear();
			for (int i = 0; i < parsedNames.size(); i++) {
				if (names.count(parsedNames[i]) == 0) {
					validSecond.push_back(parsedNames[i]);
				}
			}
			
			if ((dups) && (validSecond.size() != parsedNames.size())) {  //if dups is true and we want to get rid of anyone, get rid of everyone
				for (int i = 0; i < parsedNames.size(); i++) {  names.insert(parsedNames[i]);  }
			}else {
				//if the name in the first column is in the set then print it and any other names in second column also in set
				if (names.count(firstCol) == 0) {
					
					wroteSomething = true;
					
					out << firstCol << '\t';
					
					//you know you have at least one valid second since first column is valid
					for (int i = 0; i < validSecond.size()-1; i++) {  out << validSecond[i] << ',';  }
					out << validSecond[validSecond.size()-1] << endl;
					
					//make first name in set you come to first column and then add the remaining names to second column
				}else {
					
					//you want part of this row
					if (validSecond.size() != 0) {
						
						wroteSomething = true;
						
						out << validSecond[0] << '\t';
						
						//you know you have at least one valid second since first column is valid
						for (int i = 0; i < validSecond.size()-1; i++) {  out << validSecond[i] << ',';  }
						out << validSecond[validSecond.size()-1] << endl;
					}
				}
			}
			m->gobble(in);
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the .accnos file."); m->mothurOutEndLine();  }
		outputTypes["name"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "readName");
		exit(1);
	}
}

//**********************************************************************************************************************
int RemoveSeqsCommand::readGroup(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(groupfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(groupfile)) + "pick" + m->getExtension(groupfile);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);

		ifstream in;
		m->openInputFile(groupfile, in);
		string name, group;
		
		bool wroteSomething = false;
		
		while(!in.eof()){
			if (m->control_pressed) { in.close();  out.close();  remove(outputFileName.c_str());  return 0; }
			
			in >> name;				//read from first column
			in >> group;			//read from second column
			
			//if this name is in the accnos file
			if (names.count(name) == 0) {
				wroteSomething = true;
				out << name << '\t' << group << endl;
			}
					
			m->gobble(in);
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the .accnos file."); m->mothurOutEndLine();  }
		outputTypes["group"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "readGroup");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveSeqsCommand::readTax(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(taxfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(taxfile)) + "pick" + m->getExtension(taxfile);
		ofstream out;
		m->openOutputFile(outputFileName, out);

		ifstream in;
		m->openInputFile(taxfile, in);
		string name, tax;
		
		bool wroteSomething = false;
		
		while(!in.eof()){
			if (m->control_pressed) { in.close();  out.close();  remove(outputFileName.c_str());  return 0; }
			
			in >> name;				//read from first column
			in >> tax;			//read from second column
			
			//if this name is in the accnos file
			if (names.count(name) == 0) {
				wroteSomething = true;
				out << name << '\t' << tax << endl;
			}
					
			m->gobble(in);
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the .accnos file."); m->mothurOutEndLine();  }
		outputTypes["taxonomy"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "readTax");
		exit(1);
	}
}
//**********************************************************************************************************************
//alignreport file has a column header line then all other lines contain 16 columns.  we just want the first column since that contains the name
int RemoveSeqsCommand::readAlign(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(alignfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(alignfile)) + "pick.align.report";
		
		ofstream out;
		m->openOutputFile(outputFileName, out);

		ifstream in;
		m->openInputFile(alignfile, in);
		string name, junk;
		
		bool wroteSomething = false;
		
		//read column headers
		for (int i = 0; i < 16; i++) {  
			if (!in.eof())	{	in >> junk;	 out << junk << '\t';	}
			else			{	break;			}
		}
		out << endl;
		
		while(!in.eof()){
			if (m->control_pressed) { in.close();  out.close();  remove(outputFileName.c_str());  return 0; }
			
			in >> name;				//read from first column
			
			//if this name is in the accnos file
			if (names.count(name) == 0) {
				wroteSomething = true;
				
				out << name << '\t';
				
				//read rest
				for (int i = 0; i < 15; i++) {  
					if (!in.eof())	{	in >> junk;	 out << junk << '\t';	}
					else			{	break;			}
				}
				out << endl;
				
			}else {//still read just don't do anything with it
				
				//read rest
				for (int i = 0; i < 15; i++) {  
					if (!in.eof())	{	in >> junk;		}
					else			{	break;			}
				}
			}
			
			m->gobble(in);
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the .accnos file."); m->mothurOutEndLine();  }
		outputTypes["alignreport"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "readAlign");
		exit(1);
	}
}
//**********************************************************************************************************************
void RemoveSeqsCommand::readAccnos(){
	try {
		
		ifstream in;
		m->openInputFile(accnosfile, in);
		string name;
		
		while(!in.eof()){
			in >> name;
						
			names.insert(name);
			
			m->gobble(in);
		}
		in.close();		

	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "readAccnos");
		exit(1);
	}
}

//**********************************************************************************************************************


