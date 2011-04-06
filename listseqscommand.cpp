/*
 *  listseqscommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "listseqscommand.h"
#include "sequence.hpp"
#include "listvector.hpp"


//**********************************************************************************************************************
vector<string> ListSeqsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "FNGLT", "FNGLT", "none",false,false); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "FNGLT", "FNGLT", "none",false,false); parameters.push_back(pname);
		CommandParameter pgroup("group", "InputTypes", "", "", "FNGLT", "FNGLT", "none",false,false); parameters.push_back(pgroup);
		CommandParameter plist("list", "InputTypes", "", "", "FNGLT", "FNGLT", "none",false,false); parameters.push_back(plist);
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "FNGLT", "FNGLT", "none",false,false); parameters.push_back(ptaxonomy);
		CommandParameter palignreport("alignreport", "InputTypes", "", "", "FNGLT", "FNGLT", "none",false,false); parameters.push_back(palignreport);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ListSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ListSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The list.seqs command reads a fasta, name, group, list, taxonomy or alignreport file and outputs a .accnos file containing sequence names.\n";
		helpString += "The list.seqs command parameters are fasta, name, group, list, taxonomy and alignreport.  You must provide one of these parameters.\n";
		helpString += "The list.seqs command should be in the following format: list.seqs(fasta=yourFasta).\n";
		helpString += "Example list.seqs(fasta=amazon.fasta).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ListSeqsCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
ListSeqsCommand::ListSeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["accnos"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ListSeqsCommand", "ListSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

ListSeqsCommand::ListSeqsCommand(string option)  {
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
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["accnos"] = tempOutNames;
			
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
			}

			//check for required parameters
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
			
			if ((fastafile == "") && (namefile == "") && (listfile == "") && (groupfile == "") && (alignfile == "") && (taxfile == ""))  { m->mothurOut("You must provide a file."); m->mothurOutEndLine(); abort = true; }
			
			int okay = 1;
			if (outputDir != "") { okay++; }
			if (inputDir != "") { okay++; }
			
			if (parameters.size() > okay) { m->mothurOut("You may only enter one file."); m->mothurOutEndLine(); abort = true;  }
		}

	}
	catch(exception& e) {
		m->errorOut(e, "ListSeqsCommand", "ListSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int ListSeqsCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//read functions fill names vector
		if (fastafile != "")		{	inputFileName = fastafile;	readFasta();	}
		else if (namefile != "")	{	inputFileName = namefile;	readName();		}
		else if (groupfile != "")	{	inputFileName = groupfile;	readGroup();	}
		else if (alignfile != "")	{	inputFileName = alignfile;	readAlign();	}
		else if (listfile != "")	{	inputFileName = listfile;	readList();		}
		else if (taxfile != "")		{	inputFileName = taxfile;	readTax();		}
		
		if (m->control_pressed) { outputTypes.clear();  return 0; }
		
		//sort in alphabetical order
		sort(names.begin(), names.end());
		
		if (outputDir == "") {  outputDir += m->hasPath(inputFileName);  }
		
		string outputFileName = outputDir + m->getRootName(m->getSimpleName(inputFileName)) + "accnos";

		ofstream out;
		m->openOutputFile(outputFileName, out);
		outputNames.push_back(outputFileName); outputTypes["accnos"].push_back(outputFileName);
		
		//output to .accnos file
		for (int i = 0; i < names.size(); i++) {
			
			if (m->control_pressed) { outputTypes.clear(); out.close(); remove(outputFileName.c_str()); return 0; }
			
			out << names[i] << endl;
		}
		out.close();
		
		if (m->control_pressed) { outputTypes.clear();  remove(outputFileName.c_str()); return 0; }
		
		m->setAccnosFile(outputFileName);
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		m->mothurOut(outputFileName); m->mothurOutEndLine();	
		m->mothurOutEndLine();
		
		//set accnos file as new current accnosfile
		string current = "";
		itTypes = outputTypes.find("accnos");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setAccnosFile(current); }
		}
		
		return 0;		
	}

	catch(exception& e) {
		m->errorOut(e, "ListSeqsCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
int ListSeqsCommand::readFasta(){
	try {
		
		ifstream in;
		m->openInputFile(fastafile, in);
		string name;
		
		while(!in.eof()){
			
			if (m->control_pressed) { in.close(); return 0; }
			
			Sequence currSeq(in);
			name = currSeq.getName();
			
			if (name != "") {  names.push_back(name);  }
			
			m->gobble(in);
		}
		in.close();	
		
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "ListSeqsCommand", "readFasta");
		exit(1);
	}
}
//**********************************************************************************************************************
int ListSeqsCommand::readList(){
	try {
		ifstream in;
		m->openInputFile(listfile, in);
		
		if(!in.eof()){
			//read in list vector
			ListVector list(in);
			
			//for each bin
			for (int i = 0; i < list.getNumBins(); i++) {
				string binnames = list.get(i);
				
				if (m->control_pressed) { in.close(); return 0; }
				
				while (binnames.find_first_of(',') != -1) { 
					string name = binnames.substr(0,binnames.find_first_of(','));
					binnames = binnames.substr(binnames.find_first_of(',')+1, binnames.length());
					names.push_back(name);
				}
			
				names.push_back(binnames);
			}
		}
		in.close();	
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ListSeqsCommand", "readList");
		exit(1);
	}
}

//**********************************************************************************************************************
int ListSeqsCommand::readName(){
	try {
		
		ifstream in;
		m->openInputFile(namefile, in);
		string name, firstCol, secondCol;
		
		while(!in.eof()){
		
			if (m->control_pressed) { in.close(); return 0; }

			in >> firstCol;				
			in >> secondCol;			
			
			//parse second column saving each name
			while (secondCol.find_first_of(',') != -1) { 
				name = secondCol.substr(0,secondCol.find_first_of(','));
				secondCol = secondCol.substr(secondCol.find_first_of(',')+1, secondCol.length());
				names.push_back(name);
			}
			
			//get name after last ,
			names.push_back(secondCol);
			
			m->gobble(in);
		}
		in.close();
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ListSeqsCommand", "readName");
		exit(1);
	}
}

//**********************************************************************************************************************
int ListSeqsCommand::readGroup(){
	try {
	
		ifstream in;
		m->openInputFile(groupfile, in);
		string name, group;
		
		while(!in.eof()){
			
			if (m->control_pressed) { in.close(); return 0; }
			
			in >> name;	m->gobble(in);			//read from first column
			in >> group;			//read from second column
			
			names.push_back(name);
					
			m->gobble(in);
		}
		in.close();
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "ListSeqsCommand", "readGroup");
		exit(1);
	}
}

//**********************************************************************************************************************
//alignreport file has a column header line then all other lines contain 16 columns.  we just want the first column since that contains the name
int ListSeqsCommand::readAlign(){
	try {
	
		ifstream in;
		m->openInputFile(alignfile, in);
		string name, junk;
		
		//read column headers
		for (int i = 0; i < 16; i++) {  
			if (!in.eof())	{	in >> junk;		}
			else			{	break;			}
		}
		
		
		while(!in.eof()){
		
			if (m->control_pressed) { in.close(); return 0; }

			in >> name;				//read from first column
			
			//read rest
			for (int i = 0; i < 15; i++) {  
				if (!in.eof())	{	in >> junk;		}
				else			{	break;			}
			}
			
			names.push_back(name);
					
			m->gobble(in);
		}
		in.close();
		
		return 0;

		
	}
	catch(exception& e) {
		m->errorOut(e, "ListSeqsCommand", "readAlign");
		exit(1);
	}
}
//**********************************************************************************************************************
int ListSeqsCommand::readTax(){
	try {
		
		ifstream in;
		m->openInputFile(taxfile, in);
		string name, firstCol, secondCol;
		
		while(!in.eof()){
		
			if (m->control_pressed) { in.close(); return 0; }

			in >> firstCol;				
			in >> secondCol;			
			
			names.push_back(firstCol);
			
			m->gobble(in);
			
		}
		in.close();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ListSeqsCommand", "readTax");
		exit(1);
	}
}
//**********************************************************************************************************************
