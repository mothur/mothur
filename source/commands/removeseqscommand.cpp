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
#include "counttable.h"
#include "fastqread.h"
#include "inputdata.h"

//**********************************************************************************************************************
vector<string> RemoveSeqsCommand::setParameters(){	
	try {
        CommandParameter pfastq("fastq", "InputTypes", "", "", "none", "FNGLT", "none","fastq",false,false,true); parameters.push_back(pfastq);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "FNGLT", "none","fasta",false,false,true); parameters.push_back(pfasta);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "FNGLT", "none","name",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "FNGLT", "none","count",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "FNGLT", "none","group",false,false,true); parameters.push_back(pgroup);
		CommandParameter plist("list", "InputTypes", "", "", "none", "FNGLT", "none","list",false,false,true); parameters.push_back(plist);
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "FNGLT", "none","taxonomy",false,false,true); parameters.push_back(ptaxonomy);
		CommandParameter palignreport("alignreport", "InputTypes", "", "", "none", "FNGLT", "none","alignreport",false,false); parameters.push_back(palignreport);
		CommandParameter pqfile("qfile", "InputTypes", "", "", "none", "FNGLT", "none","qfile",false,false); parameters.push_back(pqfile);
		CommandParameter paccnos("accnos", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(paccnos);
		CommandParameter pdups("dups", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pdups);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pformat("format", "Multiple", "sanger-illumina-solexa-illumina1.8+", "illumina1.8+", "", "", "","",false,false,true); parameters.push_back(pformat);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
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
		helpString += "The remove.seqs command reads an .accnos file and at least one of the following file types: fasta, name, group, count, list, taxonomy, quality, fastq or alignreport file.\n";
		helpString += "It outputs a file containing the sequences NOT in the .accnos file.\n";
		helpString += "The remove.seqs command parameters are accnos, fasta, name, group, count, list, taxonomy, qfile, alignreport, fastq and dups.  You must provide accnos and at least one of the file parameters.\n";
        helpString += "The format parameter is used to indicate whether your sequences are sanger, solexa, illumina1.8+ or illumina, default=illumina1.8+.\n";
		helpString += "The dups parameter allows you to remove the entire line from a name file if you remove any name from the line. default=true. \n";
		helpString += "The remove.seqs command should be in the following format: remove.seqs(accnos=yourAccnos, fasta=yourFasta).\n";
		helpString += "Example remove.seqs(accnos=amazon.accnos, fasta=amazon.fasta).\n";
		;
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string RemoveSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta")            {   pattern = "[filename],pick,[extension]";    }
        else if (type == "fastq")       {   pattern = "[filename],pick,[extension]";    }
        else if (type == "taxonomy")    {   pattern = "[filename],pick,[extension]";    }
        else if (type == "name")        {   pattern = "[filename],pick,[extension]";    }
        else if (type == "group")       {   pattern = "[filename],pick,[extension]";    }
        else if (type == "count")       {   pattern = "[filename],pick,[extension]";    }
        else if (type == "list")        {   pattern = "[filename],[distance],pick,[extension]";    }
        else if (type == "qfile")       {   pattern = "[filename],pick,[extension]";    }
        else if (type == "alignreport")      {   pattern = "[filename],pick.align.report";    }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "GetSeqsCommand", "getOutputPattern");
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
        outputTypes["fastq"] = tempOutNames;
		outputTypes["taxonomy"] = tempOutNames;
		outputTypes["name"] = tempOutNames;
		outputTypes["group"] = tempOutNames;
		outputTypes["alignreport"] = tempOutNames;
		outputTypes["list"] = tempOutNames;
		outputTypes["qfile"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
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
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
            outputTypes["fastq"] = tempOutNames;
			outputTypes["taxonomy"] = tempOutNames;
			outputTypes["name"] = tempOutNames;
			outputTypes["group"] = tempOutNames;
			outputTypes["alignreport"] = tempOutNames;
			outputTypes["list"] = tempOutNames;
			outputTypes["qfile"] = tempOutNames;
            outputTypes["count"] = tempOutNames;
			
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
				
				it = parameters.find("accnos");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["accnos"] = inputDir + it->second;		}
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
				
				it = parameters.find("qfile");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["qfile"] = inputDir + it->second;		}
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
			accnosfile = validParameter.validFile(parameters, "accnos");
			if (accnosfile == "not open") { abort = true; }
			else if (accnosfile == "not found") {  
				accnosfile = current->getAccnosFile();
				if (accnosfile != "") {  m->mothurOut("Using " + accnosfile + " as input file for the accnos parameter."); m->mothurOutEndLine(); }
				else { 
					m->mothurOut("You have no valid accnos file and accnos is required."); m->mothurOutEndLine(); 
					abort = true;
				}  
			}else { current->setAccnosFile(accnosfile); }	
			
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not open") { fastafile = ""; abort = true; }
			else if (fastafile == "not found") {  fastafile = "";  }	
			else { current->setFastaFile(fastafile); }
								   
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { namefile = ""; abort = true; }
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
			
			qualfile = validParameter.validFile(parameters, "qfile");
			if (qualfile == "not open") { abort = true; }
			else if (qualfile == "not found") {  qualfile = "";  }			
			else { current->setQualFile(qualfile); }
            
            fastqfile = validParameter.validFile(parameters, "fastq");
			if (fastqfile == "not open") { abort = true; }
			else if (fastqfile == "not found") {  fastqfile = "";  }
			
			string usedDups = "true";
			string temp = validParameter.valid(parameters, "dups");
			if (temp == "not found") { 
				if (namefile != "") {  temp = "true";					}
				else				{  temp = "false"; usedDups = "";	}
			}
			dups = util.isTrue(temp);
            
            countfile = validParameter.validFile(parameters, "count");
            if (countfile == "not open") { countfile = ""; abort = true; }
            else if (countfile == "not found") { countfile = "";  }	
            else { current->setCountFile(countfile); }
            
            if ((namefile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: name or count."); m->mothurOutEndLine(); abort = true;
            }
            
            if ((groupfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: group or count."); m->mothurOutEndLine(); abort=true;
            }
			
			if ((fastqfile == "") && (countfile == "") && (fastafile == "") && (namefile == "") && (groupfile == "") && (alignfile == "") && (listfile == "") && (taxfile == "") && (qualfile == ""))  { m->mothurOut("You must provide at least one of the following: fasta, name, group, taxonomy, quality, alignreport, fastq or list."); m->mothurOutEndLine(); abort = true; }
			
            if (countfile == "") {
                if ((fastafile != "") && (namefile == "")) {
                    vector<string> files; files.push_back(fastafile);
                    if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
                }
            }
            
            format = validParameter.valid(parameters, "format");		if (format == "not found"){	format = "illumina1.8+";	}
            
            if ((format != "sanger") && (format != "illumina") && (format != "illumina1.8+") && (format != "solexa"))  {
                m->mothurOut(format + " is not a valid format. Your format choices are sanger, solexa, illumina1.8+ and illumina, aborting." ); m->mothurOutEndLine();
                abort=true;
            }
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
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		//get names you want to keep
		names = util.readAccnos(accnosfile);
		
		if (m->getControl_pressed()) { return 0; }
        
        if (countfile != "") {
            if ((fastafile != "") || (listfile != "") || (taxfile != "")) { 
                m->mothurOut("\n[NOTE]: The count file should contain only unique names, so mothur assumes your fasta, list and taxonomy files also contain only uniques.\n\n");
            }
        }
		
		//read through the correct file and output lines you want to keep
		if (namefile != "")			{		readName();		}
		if (fastafile != "")		{		readFasta();	}
        if (fastqfile != "")		{		readFastq();		}
		if (groupfile != "")		{		readGroup();	}
		if (alignfile != "")		{		readAlign();	}
		if (listfile != "")			{		readList();		}
		if (taxfile != "")			{		readTax();		}
		if (qualfile != "")			{		readQual();		}
        if (countfile != "")		{		readCount();		}
		
		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
	
		if (outputNames.size() != 0) {
			m->mothurOutEndLine();
			m->mothurOut("Output File Names: "); m->mothurOutEndLine();
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
			m->mothurOutEndLine();
			
			//set fasta file as new current fastafile
			string currentName = "";
			itTypes = outputTypes.find("fasta");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); }
			}
			
			itTypes = outputTypes.find("name");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setNameFile(currentName); }
			}
			
			itTypes = outputTypes.find("group");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setGroupFile(currentName); }
			}
			
			itTypes = outputTypes.find("list");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setListFile(currentName); }
			}
			
			itTypes = outputTypes.find("taxonomy");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setTaxonomyFile(currentName); }
			}
			
			itTypes = outputTypes.find("qfile");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setQualFile(currentName); }
			}	
            
            itTypes = outputTypes.find("count");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setCountFile(currentName); }
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
		if (outputDir == "") {  thisOutputDir += util.hasPath(fastafile);  }
		map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(fastafile));
        variables["[extension]"] = util.getExtension(fastafile);
		string outputFileName = getOutputFileName("fasta", variables);
		
		ofstream out;
		util.openOutputFile(outputFileName, out);
		
		ifstream in;
		util.openInputFile(fastafile, in);
		string name;
		
		bool wroteSomething = false;
		int removedCount = 0;
		
        set<string> uniqueNames;
		while(!in.eof()){
			if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
			
			Sequence currSeq(in);
            
            if (!dups) {//adjust name if needed
                map<string, string>::iterator it = uniqueMap.find(currSeq.getName());
                if (it != uniqueMap.end()) { currSeq.setName(it->second); }
            }

			name = currSeq.getName();
			
			if (name != "") {
				//if this name is in the accnos file
				if (names.count(name) == 0) {
                    if (uniqueNames.count(name) == 0) { //this name hasn't been seen yet
                        uniqueNames.insert(name);
                        wroteSomething = true;
					
                        currSeq.printSequence(out);
                    }else {
                        m->mothurOut("[WARNING]: " + name + " is in your fasta file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                    }
				}else {  removedCount++;  }
			}
			util.gobble(in);
		}
		in.close();	
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the .accnos file."); m->mothurOutEndLine();  }
		outputTypes["fasta"].push_back(outputFileName);  outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from your fasta file."); m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "readFasta");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveSeqsCommand::readFastq(){
	try {
		bool wroteSomething = false;
		int removedCount = 0;
        
		ifstream in;
		util.openInputFile(fastqfile, in);
		
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(fastqfile);  }
		map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(fastqfile));
        variables["[extension]"] = util.getExtension(fastqfile);
		string outputFileName = getOutputFileName("fastq", variables);
		ofstream out;
		util.openOutputFile(outputFileName, out);
        
		set<string> uniqueNames;
		while(!in.eof()){
			
			if (m->getControl_pressed()) { in.close(); out.close(); util.mothurRemove(outputFileName); return 0; }
			
            //read sequence name
            bool ignore;
            FastqRead fread(in, ignore, format); util.gobble(in);
            
            if (!ignore) {
                string name = fread.getName();
                
                if (names.count(name) == 0) {
                    if (uniqueNames.count(name) == 0) { //this name hasn't been seen yet
                        wroteSomething = true;
                        fread.printFastq(out);
                        uniqueNames.insert(name);
                    }else {
                        m->mothurOut("[WARNING]: " + name + " is in your fastq file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                    }

                }else { removedCount++; }
            }
            
			util.gobble(in);
		}
		in.close();
		out.close();
		
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the .accnos file."); m->mothurOutEndLine();  }
		outputTypes["fasta"].push_back(outputFileName);  outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from your fastq file."); m->mothurOutEndLine();

		
		return 0;
        
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "readFastq");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveSeqsCommand::readQual(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(qualfile);  }
		map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(qualfile));
        variables["[extension]"] = util.getExtension(qualfile);
		string outputFileName = getOutputFileName("qfile", variables);
		ofstream out;
		util.openOutputFile(outputFileName, out);
		
		
		ifstream in;
		util.openInputFile(qualfile, in);
		string name;
		
		bool wroteSomething = false;
		int removedCount = 0;
		
		set<string> uniqueNames;
		while(!in.eof()){	
			string saveName = "";
			string name = "";
			string scores = "";
			
			in >> name; 
			
			if (name.length() != 0) { 
				saveName = name.substr(1);
				while (!in.eof())	{	
					char c = in.get(); 
					if (c == 10 || c == 13 || c == -1){	break;	}
					else { name += c; }	
				} 
				util.gobble(in);
			}
			
			while(in){
				char letter= in.get();
				if(letter == '>'){	in.putback(letter);	break;	}
				else{ scores += letter; }
			}
			
			util.gobble(in);
			
            if (!dups) {//adjust name if needed
                map<string, string>::iterator it = uniqueMap.find(saveName);
                if (it != uniqueMap.end()) { name = ">" + it->second; saveName = it->second; }
            }
            
			if (names.count(saveName) == 0) {
                if (uniqueNames.count(saveName) == 0) { //this name hasn't been seen yet
                    uniqueNames.insert(saveName);
                    wroteSomething = true;
				
                    out << name << endl << scores;
                }else {
                    m->mothurOut("[WARNING]: " + saveName + " is in your qfile more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                }
			}else {  removedCount++;  }
			
			util.gobble(in);
		}
		in.close();
		out.close();
		
		
		if (wroteSomething == false) { m->mothurOut("Your file contains only sequences from the .accnos file."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName);  outputTypes["qfile"].push_back(outputFileName); 
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from your quality file."); m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "readQual");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveSeqsCommand::readCount(){
	try {
        
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(countfile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(countfile));
        variables["[extension]"] = util.getExtension(countfile);
		string outputFileName = getOutputFileName("count", variables);
		
        CountTable ct; ct.readTable(countfile, true, false); int originalCount = ct.getNumSeqs();
        
        for (set<string>::iterator it = names.begin(); it != names.end(); it++) { ct.remove(*it); }
        
        if (ct.getNumSeqs() == 0) {  m->mothurOut("Your file contains only sequences from the .accnos file.\n");   return 0; }
        
        ct.printTable(outputFileName);
        
        int removedCount = originalCount - ct.getNumSeqs();
        
		outputTypes["count"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from your count file.\n");
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "readCount");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveSeqsCommand::readList(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(listfile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(listfile));
        variables["[extension]"] = util.getExtension(listfile);
        InputData input(listfile, "list", nullVector);
        ListVector* list = input.getListVector();
		
		bool wroteSomething = false;
		int removedCount = 0;
        
		while(list != NULL) {
			
			removedCount = 0;
            set<string> uniqueNames;
			
            //make a new list vector
			ListVector newList;
			newList.setLabel(list->getLabel());
            
			variables["[distance]"] = list->getLabel();
            string outputFileName = getOutputFileName("list", variables);
			
			ofstream out;
			util.openOutputFile(outputFileName, out);
			outputTypes["list"].push_back(outputFileName);  outputNames.push_back(outputFileName);
            
            vector<string> binLabels = list->getLabels();
            vector<string> newBinLabels;
            
            if (m->getControl_pressed()) { out.close();  return 0; }

			//for each bin
			for (int i = 0; i < list->getNumBins(); i++) {
				if (m->getControl_pressed()) {  out.close();  util.mothurRemove(outputFileName);  return 0; }
			
				//parse out names that are in accnos file
				string bin = list->get(i);
                vector<string> bnames;
                util.splitAtComma(bin, bnames);
				
				string newNames = "";
                for (int j = 0; j < bnames.size(); j++) {
					string name = bnames[j];
                    //if that name is in the .accnos file, add it
					if (names.count(name) == 0) {
                        if (uniqueNames.count(name) == 0) { //this name hasn't been seen yet
                            uniqueNames.insert(name);
                            newNames += name + ",";
                        }else {
                            m->mothurOut("[WARNING]: " + name + " is in your list file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                        }
                    }
					else {  removedCount++;  }
                }

				//if there are names in this bin add to new list
				if (newNames != "") {  
					newNames = newNames.substr(0, newNames.length()-1); //rip off extra comma
					newList.push_back(newNames);
                    newBinLabels.push_back(binLabels[i]);
				}
			}
				
			//print new listvector
			if (newList.getNumBins() != 0) {
				wroteSomething = true;
				newList.setLabels(newBinLabels);
				newList.print(out, false);

			}

            out.close();
            
            delete list;
            list = input.getListVector();
		}
		
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the .accnos file."); m->mothurOutEndLine();  }
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from your list file."); m->mothurOutEndLine();
		
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
		if (outputDir == "") {  thisOutputDir += util.hasPath(namefile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(namefile));
        variables["[extension]"] = util.getExtension(namefile);
		string outputFileName = getOutputFileName("name", variables);
		ofstream out;
		util.openOutputFile(outputFileName, out);

		ifstream in;
		util.openInputFile(namefile, in);
		string name, firstCol, secondCol;
		
		bool wroteSomething = false;
		int removedCount = 0;
		
        set<string> uniqueNames;
		while(!in.eof()){
			if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
			
			in >> firstCol;		util.gobble(in);		
			in >> secondCol;			
			
			vector<string> parsedNames;
			util.splitAtComma(secondCol, parsedNames);
			
            vector<string> validSecond;  validSecond.clear(); vector<string> parsedNames2;
            bool parsedError = false;
			for (int i = 0; i < parsedNames.size(); i++) {
				if (names.count(parsedNames[i]) == 0) {
                    if (uniqueNames.count(parsedNames[i]) == 0) { //this name hasn't been seen yet
                        uniqueNames.insert(parsedNames[i]);
                        validSecond.push_back(parsedNames[i]);
                        parsedNames2.push_back(parsedNames[i]);
                    }else {
                        m->mothurOut("[WARNING]: " + parsedNames[i] + " is in your name file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                        parsedError = true;
                    }
				}
			}
            
            if (parsedError) {  parsedNames = parsedNames2; }
			
			if ((dups) && (validSecond.size() != parsedNames.size())) {  //if dups is true and we want to get rid of anyone, get rid of everyone
				for (int i = 0; i < parsedNames.size(); i++) {  names.insert(parsedNames[i]);  }
				removedCount += parsedNames.size();
			}else {
                if (validSecond.size() != 0) {
                    removedCount += parsedNames.size()-validSecond.size();
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
                            //we are changing the unique name in the fasta file
                            uniqueMap[firstCol] = validSecond[0];
                            
                            //you know you have at least one valid second since first column is valid
                            for (int i = 0; i < validSecond.size()-1; i++) {  out << validSecond[i] << ',';  }
                            out << validSecond[validSecond.size()-1] << endl;
                        }
                    }
                }
			}
			util.gobble(in);
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the .accnos file."); m->mothurOutEndLine();  }
		outputTypes["name"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from your name file."); m->mothurOutEndLine();
		
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
		if (outputDir == "") {  thisOutputDir += util.hasPath(groupfile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(groupfile));
        variables["[extension]"] = util.getExtension(groupfile);
		string outputFileName = getOutputFileName("group", variables);	
		ofstream out;
		util.openOutputFile(outputFileName, out);

		ifstream in;
		util.openInputFile(groupfile, in);
		string name, group;
		
		bool wroteSomething = false;
		int removedCount = 0;
		
        set<string> uniqueNames;
		while(!in.eof()){
			if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
			
			in >> name;			util.gobble(in);		//read from first column
			in >> group;			//read from second column
			
			//if this name is in the accnos file
			if (names.count(name) == 0) {
                if (uniqueNames.count(name) == 0) { //this name hasn't been seen yet
                    uniqueNames.insert(name);
                    wroteSomething = true;
                    out << name << '\t' << group << endl;
                }else {
                    m->mothurOut("[WARNING]: " + name + " is in your group file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                }
			}else {  removedCount++;  }
					
			util.gobble(in);
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the .accnos file."); m->mothurOutEndLine();  }
		outputTypes["group"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from your group file."); m->mothurOutEndLine();

		
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
		if (outputDir == "") {  thisOutputDir += util.hasPath(taxfile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(taxfile));
        variables["[extension]"] = util.getExtension(taxfile);
		string outputFileName = getOutputFileName("taxonomy", variables);
		ofstream out;
		util.openOutputFile(outputFileName, out);

		ifstream in;
		util.openInputFile(taxfile, in);
		string name, tax;
		
		bool wroteSomething = false;
		int removedCount = 0;
		
        set<string> uniqueNames;
		while(!in.eof()){
			if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
			
            in >> name; util.gobble(in);
            tax = util.getline(in); util.gobble(in);
            
            if (!dups) {//adjust name if needed
                map<string, string>::iterator it = uniqueMap.find(name);
                if (it != uniqueMap.end()) { name = it->second; }
            }
            
			//if this name is in the accnos file
			if (names.count(name) == 0) {
                if (uniqueNames.count(name) == 0) { //this name hasn't been seen yet
                    uniqueNames.insert(name);
                    wroteSomething = true;
            
                    out << name << '\t' << tax << endl;
                }else {
                    m->mothurOut("[WARNING]: " + name + " is in your taxonomy file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                }
			}else {  removedCount++;  }
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the .accnos file."); m->mothurOutEndLine();  }
		outputTypes["taxonomy"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from your taxonomy file."); m->mothurOutEndLine();
		
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
		if (outputDir == "") {  thisOutputDir += util.hasPath(alignfile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(alignfile));
		string outputFileName = getOutputFileName("alignreport", variables);
		
		ofstream out;
		util.openOutputFile(outputFileName, out);

		ifstream in;
		util.openInputFile(alignfile, in);
		string name, junk;
		
		bool wroteSomething = false;
		int removedCount = 0;
		
		//read column headers
		for (int i = 0; i < 16; i++) {  
			if (!in.eof())	{	in >> junk;	 out << junk << '\t';	}
			else			{	break;			}
		}
		out << endl;
		
        set<string> uniqueNames;
		while(!in.eof()){
			if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
			
			in >> name;				//read from first column
            
            if (!dups) {//adjust name if needed
                map<string, string>::iterator it = uniqueMap.find(name);
                if (it != uniqueMap.end()) { name = it->second; }
            }
			
			//if this name is in the accnos file
			if (names.count(name) == 0) {
                if (uniqueNames.count(name) == 0) { //this name hasn't been seen yet
                    uniqueNames.insert(name);
                    wroteSomething = true;
                    
                    out << name << '\t';
                    
                    //read rest
                    for (int i = 0; i < 15; i++) {
                        if (!in.eof())	{	in >> junk;	 out << junk << '\t';	}
                        else			{	break;			}
                    }
                    out << endl;
                }else {
                    m->mothurOut("[WARNING]: " + name + " is in your alignreport file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                }
			}else {//still read just don't do anything with it
				removedCount++;  
				
				//read rest
				for (int i = 0; i < 15; i++) {  
					if (!in.eof())	{	in >> junk;		}
					else			{	break;			}
				}
			}
			
			util.gobble(in);
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the .accnos file."); m->mothurOutEndLine();  }
		outputTypes["alignreport"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from your alignreport file."); m->mothurOutEndLine();

		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveSeqsCommand", "readAlign");
		exit(1);
	}
}
//**********************************************************************************************************************


