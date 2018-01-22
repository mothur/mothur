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
#include "counttable.h"
#include "fastqread.h"

//**********************************************************************************************************************
vector<string> ListSeqsCommand::setParameters(){	
	try {
        CommandParameter pfastq("fastq", "InputTypes", "", "", "FNGLT", "FNGLT", "none","accnos",false,false,true); parameters.push_back(pfastq);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "FNGLT", "FNGLT", "none","accnos",false,false,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "FNGLT", "FNGLT", "none","accnos",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "FNGLT", "FNGLT", "none","accnos",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "FNGLT", "FNGLT", "none","accnos",false,false,true); parameters.push_back(pgroup);
		CommandParameter plist("list", "InputTypes", "", "", "FNGLT", "FNGLT", "none","accnos",false,false,true); parameters.push_back(plist);
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "FNGLT", "FNGLT", "none","accnos",false,false,true); parameters.push_back(ptaxonomy);
		CommandParameter palignreport("alignreport", "InputTypes", "", "", "FNGLT", "FNGLT", "none","accnos",false,false); parameters.push_back(palignreport);
        CommandParameter pformat("format", "Multiple", "sanger-illumina-solexa-illumina1.8+", "illumina1.8+", "", "", "","",false,false,true); parameters.push_back(pformat);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
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
		helpString += "The list.seqs command reads a fasta, name, group, count, list, taxonomy, fastq or alignreport file and outputs a .accnos file containing sequence names.\n";
		helpString += "The list.seqs command parameters are fasta, name, group, count, list, taxonomy, fastq and alignreport.  You must provide one of these parameters.\n";
         helpString += "The format parameter is used to indicate whether your sequences are sanger, solexa, illumina1.8+ or illumina, default=illumina1.8+.\n";
		helpString += "The list.seqs command should be in the following format: list.seqs(fasta=yourFasta).\n";
		helpString += "Example list.seqs(fasta=amazon.fasta).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ListSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ListSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "accnos") {  pattern = "[filename],accnos"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ListSeqsCommand", "getOutputPattern");
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
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
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
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";		}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("alignreport");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["alignreport"] = inputDir + it->second;		}
				}
				
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
				
				it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
                
                it = parameters.find("fastq");
				//user has given a template file
				if(it != parameters.end()){
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fastq"] = inputDir + it->second;		}
				}
			}

			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") {  fastafile = "";  }
			else { current->setFastaFile(fastafile); }
			
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { abort = true; }
			else if (namefile == "not found") {  namefile = "";  }	
			else { current->setNameFile(namefile); }
			
			groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not open") { abort = true; }
			else if (groupfile == "not found") {  groupfile = "";  }	
			else { current->setGroupFile(groupfile); }
			
			alignfile = validParameter.validFile(parameters, "alignreport");
			if (alignfile == "not open") { abort = true; }
			else if (alignfile == "not found") {  alignfile = "";  }
			
			listfile = validParameter.validFile(parameters, "list");
			if (listfile == "not open") { abort = true; }
			else if (listfile == "not found") {  listfile = "";  }
			else { current->setListFile(listfile); }
			
			taxfile = validParameter.validFile(parameters, "taxonomy");
			if (taxfile == "not open") { abort = true; }
			else if (taxfile == "not found") {  taxfile = "";  }
			else { current->setTaxonomyFile(taxfile); }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { abort = true; }
			else if (countfile == "not found") {  countfile = "";  }
			else { current->setCountFile(countfile); }
            
            fastqfile = validParameter.validFile(parameters, "fastq");
			if (fastqfile == "not open") { abort = true; }
			else if (fastqfile == "not found") {  fastqfile = "";  }
			
			if ((fastqfile == "") && (countfile == "") && (fastafile == "") && (namefile == "") && (listfile == "") && (groupfile == "") && (alignfile == "") && (taxfile == ""))  { m->mothurOut("You must provide a file."); m->mothurOutEndLine(); abort = true; }
            
            bool formatFound = true;
            format = validParameter.valid(parameters, "format");		if (format == "not found"){	formatFound = false; format = "illumina1.8+";	}
            
            if ((format != "sanger") && (format != "illumina") && (format != "illumina1.8+") && (format != "solexa"))  {
                m->mothurOut(format + " is not a valid format. Your format choices are sanger, solexa, illumina1.8+ and illumina, aborting." ); m->mothurOutEndLine();
                abort=true;
            }
            
            int okay = 1;
            if (outputDir != "") { okay++; }
            if (inputDir != "") { okay++; }
            if (formatFound) { okay++; }
			
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
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		//read functions fill names vector
		if (fastafile != "")		{	inputFileName = fastafile;	readFasta();	}
        else if (fastqfile != "")	{	inputFileName = fastqfile;	readFastq();	}
		else if (namefile != "")	{	inputFileName = namefile;	readName();		}
		else if (groupfile != "")	{	inputFileName = groupfile;	readGroup();	}
		else if (alignfile != "")	{	inputFileName = alignfile;	readAlign();	}
		else if (listfile != "")	{	inputFileName = listfile;	readList();		}
		else if (taxfile != "")		{	inputFileName = taxfile;	readTax();		}
        else if (countfile != "")	{	inputFileName = countfile;	readCount();	}
		
		if (m->getControl_pressed()) { outputTypes.clear();  return 0; }
		
		//sort in alphabetical order
		sort(names.begin(), names.end());
		
		if (outputDir == "") {  outputDir += util.hasPath(inputFileName);  }
		
        map<string, string> variables; 
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(inputFileName));
		string outputFileName = getOutputFileName("accnos", variables);

		ofstream out;
		util.openOutputFile(outputFileName, out);
		outputNames.push_back(outputFileName); outputTypes["accnos"].push_back(outputFileName);
		
		//output to .accnos file
		for (int i = 0; i < names.size(); i++) {
			
			if (m->getControl_pressed()) { outputTypes.clear(); out.close(); util.mothurRemove(outputFileName); return 0; }
			
			out << names[i] << endl;
		}
		out.close();
		
		if (m->getControl_pressed()) { outputTypes.clear();  util.mothurRemove(outputFileName); return 0; }
		
		current->setAccnosFile(outputFileName);
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		m->mothurOut(outputFileName); m->mothurOutEndLine();	
		m->mothurOutEndLine();
		
		//set accnos file as new current accnosfile
		string currentName = "";
		itTypes = outputTypes.find("accnos");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setAccnosFile(currentName); }
		}
		
		return 0;		
	}

	catch(exception& e) {
		m->errorOut(e, "ListSeqsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int ListSeqsCommand::readFastq(){
	try {
		
		ifstream in;
		util.openInputFile(fastqfile, in);
		string name;
		
		int count = 1;
		while(!in.eof()){
			
			if (m->getControl_pressed()) { in.close(); return 0; }
			
            bool ignore;
            FastqRead fread(in, ignore, format); util.gobble(in);
            
            if (!ignore) { names.push_back(fread.getName()); }
			
			if (m->getDebug()) { count++; m->mothurOut("[DEBUG]: count = " + toString(count) + ", name = " + name + "\n"); }
		}
		in.close();
		
		return 0;
        
	}
	catch(exception& e) {
		m->errorOut(e, "ListSeqsCommand", "readFastq");
		exit(1);
	}
}

//**********************************************************************************************************************
int ListSeqsCommand::readFasta(){
	try {
		
		ifstream in;
		util.openInputFile(fastafile, in);
		string name;
		
		//ofstream out;
		//string newFastaName = outputDir + util.getRootName(util.getSimpleName(fastafile)) + "numsAdded.fasta";
		//util.openOutputFile(newFastaName, out);
		int count = 1;
		//string lastName = "";
		
		while(!in.eof()){
			
			if (m->getControl_pressed()) { in.close(); return 0; }
			
			Sequence currSeq(in);
			name = currSeq.getName();
			
			if (name != "") {  names.push_back(name);  }
			
			util.gobble(in);
			if (m->getDebug()) { count++; cout << "[DEBUG]: count = " + toString(count) + ", name = " + currSeq.getName() + "\n"; }
		}
		in.close();	
		//out.close();
		
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
		util.openInputFile(listfile, in);
        string otuTag = util.getTag(listfile); //looks at filename to determine if OTU labels should be "Otu" or "Phylotype"
        string readHeaders = ""; //Tells mothur to try and read headers from the file
		
		if(!in.eof()){
			//read in list vector
			ListVector list(in, readHeaders, otuTag);
			
			//for each bin
			for (int i = 0; i < list.getNumBins(); i++) {
				string binnames = list.get(i);
				
				if (m->getControl_pressed()) { in.close(); return 0; }
				
				util.splitAtComma(binnames, names);
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
		util.openInputFile(namefile, in);
		string name, firstCol, secondCol;
		
		while(!in.eof()){
		
			if (m->getControl_pressed()) { in.close(); return 0; }

			in >> firstCol;	util.gobble(in);
			in >> secondCol;			
			
			//parse second column saving each name
			util.splitAtComma(secondCol, names);
			
			util.gobble(in);
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
		util.openInputFile(groupfile, in);
		string name, group;
		
		while(!in.eof()){
			
			if (m->getControl_pressed()) { in.close(); return 0; }
			
			in >> name;	util.gobble(in);			//read from first column
			in >> group;			//read from second column
			
			names.push_back(name);
					
			util.gobble(in);
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
int ListSeqsCommand::readCount(){
	try {
		CountTable ct;
		ct.readTable(countfile, false, false);
        
        if (m->getControl_pressed()) { return 0; }
        
        names = ct.getNamesOfSeqs();
        
        return 0;
        
	}
	catch(exception& e) {
		m->errorOut(e, "ListSeqsCommand", "readCount");
		exit(1);
	}
}
//**********************************************************************************************************************
//alignreport file has a column header line then all other lines contain 16 columns.  we just want the first column since that contains the name
int ListSeqsCommand::readAlign(){
	try {
	
		ifstream in;
		util.openInputFile(alignfile, in);
		string name, junk;
		
		//read column headers
		for (int i = 0; i < 16; i++) {  
			if (!in.eof())	{	in >> junk;		}
			else			{	break;			}
		}
		//util.getline(in);
		
		while(!in.eof()){
		
			if (m->getControl_pressed()) { in.close(); return 0; }

			in >> name;				//read from first column
			//util.getline(in);
			//read rest
			for (int i = 0; i < 15; i++) {  
				if (!in.eof())	{	in >> junk;		}
				else			{	break;			}
			}
			
			names.push_back(name);
					
			util.gobble(in);
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
		util.openInputFile(taxfile, in);
		string name, firstCol, secondCol;
		
		while(!in.eof()){
		
			if (m->getControl_pressed()) { in.close(); return 0; }

            in >> firstCol; util.gobble(in);
            secondCol = util.getline(in); util.gobble(in);
			
			names.push_back(firstCol);
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
