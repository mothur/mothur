/*
 *  getseqscommand.cpp
 *  Mothur
 *
 *  Created by Sarah Westcott on 7/8/09.
 *  Copyright 2009 Schloss Lab UMASS Amherst. All rights reserved.
 *
 */

#include "getseqscommand.h"
#include "sequence.hpp"
#include "listvector.hpp"
#include "counttable.h"

//**********************************************************************************************************************
vector<string> GetSeqsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "FNGLT", "none","fasta",false,false,true); parameters.push_back(pfasta);
        CommandParameter pfastq("fastq", "InputTypes", "", "", "none", "FNGLT", "none","fastq",false,false,true); parameters.push_back(pfastq);
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
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		CommandParameter paccnos2("accnos2", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(paccnos2);

		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The get.seqs command reads an .accnos file and any of the following file types: fasta, name, group, count, list, taxonomy, quality, fastq or alignreport file.\n";
		helpString += "It outputs a file containing only the sequences in the .accnos file.\n";
		helpString += "The get.seqs command parameters are accnos, fasta, name, group, list, taxonomy, qfile, alignreport, fastq and dups.  You must provide accnos unless you have a valid current accnos file, and at least one of the other parameters.\n";
		helpString += "The dups parameter allows you to add the entire line from a name file if you add any name from the line. default=true. \n";
		helpString += "The get.seqs command should be in the following format: get.seqs(accnos=yourAccnos, fasta=yourFasta).\n";
		helpString += "Example get.seqs(accnos=amazon.accnos, fasta=amazon.fasta).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "GetSeqsCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
GetSeqsCommand::GetSeqsCommand(){	
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
		outputTypes["accnosreport"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "GetSeqsCommand", "GetSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetSeqsCommand::getOutputPattern(string type) {
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
        else if (type == "accnosreport")      {   pattern = "[filename],pick.accnos.report";    }
        else if (type == "alignreport")      {   pattern = "[filename],pick.align.report";    }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "GetSeqsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
GetSeqsCommand::GetSeqsCommand(string option)  {
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
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
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
			outputTypes["accnosreport"] = tempOutNames;
            outputTypes["count"] = tempOutNames;
			
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
				
				it = parameters.find("accnos2");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["accnos2"] = inputDir + it->second;		}
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
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
                
                it = parameters.find("fastq");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fastq"] = inputDir + it->second;		}
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
			}else { m->setAccnosFile(accnosfile); }	
			
			if (accnosfile2 == "not found") { accnosfile2 = ""; }
			
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { fastafile = ""; abort = true; }
			else if (fastafile == "not found") {  fastafile = "";  }
			else { m->setFastaFile(fastafile); }
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") {  namefile = "";  }	
			else { m->setNameFile(namefile); }
			
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }
			else if (groupfile == "not found") {  groupfile = "";  }	
			else { m->setGroupFile(groupfile); }
			
			alignfile = validParameter.validFile(parameters, "alignreport", true);
			if (alignfile == "not open") { abort = true; }
			else if (alignfile == "not found") {  alignfile = "";  }
			
			listfile = validParameter.validFile(parameters, "list", true);
			if (listfile == "not open") { abort = true; }
			else if (listfile == "not found") {  listfile = "";  }
			else { m->setListFile(listfile); }
			
			taxfile = validParameter.validFile(parameters, "taxonomy", true);
			if (taxfile == "not open") { taxfile = ""; abort = true; }
			else if (taxfile == "not found") {  taxfile = "";  }
			else { m->setTaxonomyFile(taxfile); }
			
			qualfile = validParameter.validFile(parameters, "qfile", true);
			if (qualfile == "not open") { abort = true; }
			else if (qualfile == "not found") {  qualfile = "";  }
			else { m->setQualFile(qualfile); }
            
            fastqfile = validParameter.validFile(parameters, "fastq", true);
			if (fastqfile == "not open") { abort = true; }
			else if (fastqfile == "not found") {  fastqfile = "";  }
			
			accnosfile2 = validParameter.validFile(parameters, "accnos2", true);
			if (accnosfile2 == "not open") { abort = true; }
			else if (accnosfile2 == "not found") {  accnosfile2 = "";  }
			
            countfile = validParameter.validFile(parameters, "count", true);
            if (countfile == "not open") { countfile = ""; abort = true; }
            else if (countfile == "not found") { countfile = "";  }	
            else { m->setCountTableFile(countfile); }
            
            if ((namefile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: name or count."); m->mothurOutEndLine(); abort = true;
            }
            
            if ((groupfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: group or count."); m->mothurOutEndLine(); abort=true;
            }

			
			string usedDups = "true";
			string temp = validParameter.validFile(parameters, "dups", false);	if (temp == "not found") { temp = "true"; usedDups = ""; }
			dups = m->isTrue(temp);
			
			if ((fastqfile == "") && (fastafile == "") && (namefile == "") && (groupfile == "") && (alignfile == "") && (listfile == "") && (taxfile == "") && (qualfile == "") && (accnosfile2 == "") && (countfile == ""))  { m->mothurOut("You must provide one of the following: fasta, name, group, count, alignreport, taxonomy, quality, fastq or listfile."); m->mothurOutEndLine(); abort = true; }
            
            if (countfile == "") {
                if ((namefile == "") && ((fastafile != "") || (taxfile != ""))){
                    vector<string> files; files.push_back(fastafile); files.push_back(taxfile);
                    parser.getNameFile(files);
                }
            }
		}

	}
	catch(exception& e) {
		m->errorOut(e, "GetSeqsCommand", "GetSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetSeqsCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//get names you want to keep
		names = m->readAccnos(accnosfile);
		
		if (m->control_pressed) { return 0; }
        
        if (countfile != "") {
            if ((fastafile != "") || (listfile != "") || (taxfile != "")) { 
                m->mothurOut("\n[NOTE]: The count file should contain only unique names, so mothur assumes your fasta, list and taxonomy files also contain only uniques.\n\n");
            }
        }
		
		//read through the correct file and output lines you want to keep
		if (namefile != "")			{		readName();			}
		if (fastafile != "")		{		readFasta();		}
        if (fastqfile != "")		{		readFastq();		}
		if (groupfile != "")		{		readGroup();		}
        if (countfile != "")		{		readCount();		}
		if (alignfile != "")		{		readAlign();		}
		if (listfile != "")			{		readList();			}
		if (taxfile != "")			{		readTax();			}
		if (qualfile != "")			{		readQual();			}
		if (accnosfile2 != "")		{		compareAccnos();	}
        
        if (m->debug) { runSanityCheck(); }
		
		if (m->control_pressed) { outputTypes.clear(); for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]);  } return 0; }
		
		
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
			
            itTypes = outputTypes.find("count");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setCountTableFile(current); }
			}
		}
		
		return 0;		
	}

	catch(exception& e) {
		m->errorOut(e, "GetSeqsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int GetSeqsCommand::readFastq(){
	try {
		bool wroteSomething = false;
		int selectedCount = 0;
        
		ifstream in;
		m->openInputFile(fastqfile, in);
		
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(fastqfile);  }
		map<string, string> variables;
        variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(fastqfile));
        variables["[extension]"] = m->getExtension(fastqfile);
		string outputFileName = getOutputFileName("fastq", variables);
		ofstream out;
		m->openOutputFile(outputFileName, out);

        set<string> uniqueNames;
		while(!in.eof()){
			
			if (m->control_pressed) { in.close(); out.close(); m->mothurRemove(outputFileName); return 0; }
			
			//read sequence name
			string input = m->getline(in); m->gobble(in);
			
            string outputString = input + "\n";
            
			if (input[0] == '@') {
                //get rest of lines
                outputString += m->getline(in) + "\n"; m->gobble(in);
                outputString += m->getline(in) + "\n"; m->gobble(in);
                outputString += m->getline(in) + "\n"; m->gobble(in);
                
                vector<string> splits = m->splitWhiteSpace(input);
                string name = splits[0];
                name = name.substr(1);
                m->checkName(name);
                
                if (names.count(name) != 0) {
                    if (uniqueNames.count(name) == 0) { //this name hasn't been seen yet
                        wroteSomething = true;
                        selectedCount++;
                        out << outputString;
                        uniqueNames.insert(name);
                    }else {
                        m->mothurOut("[WARNING]: " + name + " is in your fastq file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                    }
                }
            }
            
			m->gobble(in);
		}
		in.close();
		out.close();
        
		if (wroteSomething == false) { m->mothurOut("Your file does not contain any sequence from the .accnos file."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName);  outputTypes["fastq"].push_back(outputFileName);
		
		m->mothurOut("Selected " + toString(selectedCount) + " sequences from your fastq file."); m->mothurOutEndLine();
		
		return 0;
        
	}
	catch(exception& e) {
		m->errorOut(e, "GetSeqsCommand", "readFastq");
		exit(1);
	}
}

//**********************************************************************************************************************
int GetSeqsCommand::readFasta(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(fastafile);  }
		map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(fastafile));
        variables["[extension]"] = m->getExtension(fastafile);
		string outputFileName = getOutputFileName("fasta", variables);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(fastafile, in);
		string name;
		
		bool wroteSomething = false;
		int selectedCount = 0;
        
        if (m->debug) { set<string> temp; sanity["fasta"] = temp; }
		
        set<string> uniqueNames;
		while(!in.eof()){
		
			if (m->control_pressed) { in.close(); out.close(); m->mothurRemove(outputFileName);  return 0; }
			
			Sequence currSeq(in);
			name = currSeq.getName();
            
            if (!dups) {//adjust name if needed
                map<string, string>::iterator it = uniqueMap.find(name);
                if (it != uniqueMap.end()) { currSeq.setName(it->second); }
            }
			
            name = currSeq.getName();
            
			if (name != "") {
				//if this name is in the accnos file
				if (names.count(name) != 0) {
                    if (uniqueNames.count(name) == 0) { //this name hasn't been seen yet
                        wroteSomething = true;
					
                        currSeq.printSequence(out);
                        selectedCount++;
                        uniqueNames.insert(name);
                    
                        if (m->debug) { sanity["fasta"].insert(name); }
                    }else {
                        m->mothurOut("[WARNING]: " + name + " is in your fasta file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                    }
                }
			}
			m->gobble(in);
		}
		in.close();	
		out.close();
		
		
		if (wroteSomething == false) { m->mothurOut("Your file does not contain any sequence from the .accnos file."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName);  outputTypes["fasta"].push_back(outputFileName); 
		
		m->mothurOut("Selected " + toString(selectedCount) + " sequences from your fasta file."); m->mothurOutEndLine();
		
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "GetSeqsCommand", "readFasta");
		exit(1);
	}
}
//**********************************************************************************************************************
int GetSeqsCommand::readQual(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(qualfile);  }
		map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(qualfile));
        variables["[extension]"] = m->getExtension(qualfile);
		string outputFileName = getOutputFileName("qfile", variables);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		
		ifstream in;
		m->openInputFile(qualfile, in);
		string name;
		
		bool wroteSomething = false;
		int selectedCount = 0;
		
        if (m->debug) { set<string> temp; sanity["qual"] = temp; }
		
        set<string> uniqueNames;
		while(!in.eof()){	
			string saveName = "";
			string name = "";
			string scores = "";
			
			in >> name;
            
            if (!dups) {//adjust name if needed
                map<string, string>::iterator it = uniqueMap.find(name);
                if (it != uniqueMap.end()) { name = it->second; }
            }
				
			if (name.length() != 0) { 
				saveName = name.substr(1);
				while (!in.eof())	{	
					char c = in.get(); 
					if (c == 10 || c == 13 || c == -1){	break;	}
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
			
            if (names.count(saveName) != 0) {
                if (uniqueNames.count(saveName) == 0) { //this name hasn't been seen yet
                    uniqueNames.insert(saveName);
                    wroteSomething = true;
                    
                    out << name << endl << scores;
                    selectedCount++;
                    if (m->debug) { sanity["qual"].insert(name); }
                }else {
                    m->mothurOut("[WARNING]: " + saveName + " is in your qfile more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                }
            }
			
			m->gobble(in);
		}
		in.close();
		out.close();
		
		
		if (wroteSomething == false) { m->mothurOut("Your file does not contain any sequence from the .accnos file."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName);  outputTypes["qfile"].push_back(outputFileName); 
		
		m->mothurOut("Selected " + toString(selectedCount) + " sequences from your quality file."); m->mothurOutEndLine();

		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "GetSeqsCommand", "readQual");
		exit(1);
	}
}
//**********************************************************************************************************************
int GetSeqsCommand::readCount(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(countfile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(countfile));
        variables["[extension]"] = m->getExtension(countfile);
		string outputFileName = getOutputFileName("count", variables);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(countfile, in);
		
		bool wroteSomething = false;
		int selectedCount = 0;
		
        string headers = m->getline(in); m->gobble(in);
        out << headers << endl;
        string test = headers; vector<string> pieces = m->splitWhiteSpace(test);
        
        string name, rest; int thisTotal; rest = "";
        set<string> uniqueNames;
        while (!in.eof()) {
            
            if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
            
            in >> name; m->gobble(in); 
            in >> thisTotal; m->gobble(in);
            if (pieces.size() > 2) {  rest = m->getline(in); m->gobble(in);  }
            if (m->debug) { m->mothurOut("[DEBUG]: " + name + '\t' + rest + "\n"); }
            
            if (names.count(name) != 0) {
                if (uniqueNames.count(name) == 0) { //this name hasn't been seen yet
                    uniqueNames.insert(name);

                    out << name << '\t' << thisTotal << '\t' << rest << endl;
                    wroteSomething = true;
                    selectedCount+= thisTotal;
                }else {
                    m->mothurOut("[WARNING]: " + name + " is in your count file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                }
            }
        }
        in.close();
		out.close();
        
        //check for groups that have been eliminated
        CountTable ct;
        if (ct.testGroups(outputFileName)) {
            ct.readTable(outputFileName, true, false);
            ct.printTable(outputFileName);
        }
		
		if (wroteSomething == false) {  m->mothurOut("Your file does not contain any sequence from the .accnos file."); m->mothurOutEndLine();  }
		outputTypes["count"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		m->mothurOut("Selected " + toString(selectedCount) + " sequences from your count file."); m->mothurOutEndLine();
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "GetSeqsCommand", "readCount");
		exit(1);
	}
}

//**********************************************************************************************************************
int GetSeqsCommand::readList(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(listfile);  }
        map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(listfile));
        variables["[extension]"] = m->getExtension(listfile);
		
		ifstream in;
		m->openInputFile(listfile, in);
		
		bool wroteSomething = false;
		int selectedCount = 0;
        
        if (m->debug) { set<string> temp; sanity["list"] = temp; }
		
        set<string> uniqueNames;
		while(!in.eof()){
			
			selectedCount = 0;

			//read in list vector
			ListVector list(in);
			
			//make a new list vector
			ListVector newList;
			newList.setLabel(list.getLabel());
            
            variables["[distance]"] = list.getLabel();
            string outputFileName = getOutputFileName("list", variables);
			
			ofstream out;
			m->openOutputFile(outputFileName, out);
			outputTypes["list"].push_back(outputFileName);  outputNames.push_back(outputFileName);
            
            vector<string> binLabels = list.getLabels();
            vector<string> newBinLabels;
            
            if (m->control_pressed) { in.close(); out.close();  return 0; }
			
			//for each bin
			for (int i = 0; i < list.getNumBins(); i++) {
			
				//parse out names that are in accnos file
				string binnames = list.get(i);
                vector<string> bnames;
                m->splitAtComma(binnames, bnames);
				
				string newNames = "";
                for (int j = 0; j < bnames.size(); j++) {
                    string name = bnames[j];
                    
                    //if that name is in the .accnos file, add it
                    if (names.count(name) != 0) {
                        if (uniqueNames.count(name) == 0) { //this name hasn't been seen yet
                            uniqueNames.insert(name);
                            newNames += name + ",";
                            selectedCount++;
                            if (m->debug) { sanity["list"].insert(name); }
                        }else {
                            m->mothurOut("[WARNING]: " + name + " is in your list file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                        }
                    }
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
                newList.printHeaders(out);
				newList.print(out);
			}
			
			m->gobble(in);
            out.close();
		}
		in.close();	
		
		
		if (wroteSomething == false) { m->mothurOut("Your file does not contain any sequence from the .accnos file."); m->mothurOutEndLine();  }
		
		m->mothurOut("Selected " + toString(selectedCount) + " sequences from your list file."); m->mothurOutEndLine();
		
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "GetSeqsCommand", "readList");
		exit(1);
	}
}
//**********************************************************************************************************************
int GetSeqsCommand::readName(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(namefile);  }
        map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(namefile));
        variables["[extension]"] = m->getExtension(namefile);
		string outputFileName = getOutputFileName("name", variables);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		

		ifstream in;
		m->openInputFile(namefile, in);
		string name, firstCol, secondCol;
		
		bool wroteSomething = false;
		int selectedCount = 0;
        
        if (m->debug) { set<string> temp; sanity["name"] = temp; }
        if (m->debug) { set<string> temp; sanity["dupname"] = temp; }
		
        set<string> uniqueNames;
		while(!in.eof()){
		
			if (m->control_pressed) { in.close(); out.close(); m->mothurRemove(outputFileName);  return 0; }

			in >> firstCol;			m->gobble(in);
			in >> secondCol;
			
			string hold = "";
			if (dups) { hold = secondCol; }
			
			vector<string> parsedNames;
			m->splitAtComma(secondCol, parsedNames);
			
			vector<string> validSecond; vector<string> parsedNames2;
            bool parsedError = false;
			for (int i = 0; i < parsedNames.size(); i++) {
                if (names.count(parsedNames[i]) != 0) {
                    if (uniqueNames.count(parsedNames[i]) == 0) { //this name hasn't been seen yet
                        uniqueNames.insert(parsedNames[i]);
                        validSecond.push_back(parsedNames[i]);
                        parsedNames2.push_back(parsedNames[i]);
                        if (m->debug) { sanity["dupname"].insert(parsedNames[i]); }
                    }else {
                        m->mothurOut("[WARNING]: " + parsedNames[i] + " is in your name file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                        parsedError = true;
                    }
				}
			}
            if (parsedError) {
                parsedNames = parsedNames2;
                hold = "";
                if (parsedNames.size() != 0) {
                    for (int i = 0; i < parsedNames.size()-1; i++) {  hold += parsedNames[i] + ',';  }
                    hold += parsedNames[parsedNames.size()-1] + '\n';
                }
            }

			if ((dups) && (validSecond.size() != 0)) { //dups = true and we want to add someone, then add everyone
				for (int i = 0; i < parsedNames.size(); i++) {  names.insert(parsedNames[i]); if (m->debug) { sanity["dupname"].insert(parsedNames[i]); } }
				out << firstCol << '\t' << hold << endl;
				wroteSomething = true;
				selectedCount += parsedNames.size();
                if (m->debug) { sanity["name"].insert(firstCol); }
			}else {
                
                if (validSecond.size() != 0) {
                    selectedCount += validSecond.size();
                    
                    //if the name in the first column is in the set then print it and any other names in second column also in set
                    if (names.count(firstCol) != 0) {
                        
                        wroteSomething = true;
                        
                        out << firstCol << '\t';
                        
                        //you know you have at least one valid second since first column is valid
                        for (int i = 0; i < validSecond.size()-1; i++) {  out << validSecond[i] << ',';  }
                        out << validSecond[validSecond.size()-1] << endl;
                        
                        if (m->debug) { sanity["name"].insert(firstCol); }
                        
                        
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
                            
                            if (m->debug) { sanity["name"].insert(validSecond[0]); }
                        }
                    }
                }
			}
			m->gobble(in);
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) { m->mothurOut("Your file does not contain any sequence from the .accnos file."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName); outputTypes["name"].push_back(outputFileName);
		
		m->mothurOut("Selected " + toString(selectedCount) + " sequences from your name file."); m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "GetSeqsCommand", "readName");
		exit(1);
	}
}

//**********************************************************************************************************************
int GetSeqsCommand::readGroup(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(groupfile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(groupfile));
        variables["[extension]"] = m->getExtension(groupfile);
		string outputFileName = getOutputFileName("group", variables);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		

		ifstream in;
		m->openInputFile(groupfile, in);
		string name, group;
		
		bool wroteSomething = false;
		int selectedCount = 0;
        
        if (m->debug) { set<string> temp; sanity["group"] = temp; }
		
        set<string> uniqueNames;
		while(!in.eof()){

			if (m->control_pressed) { in.close(); out.close(); m->mothurRemove(outputFileName);  return 0; }


			in >> name;				//read from first column
			in >> group;			//read from second column
            
			
            if (names.count(name) != 0) {
                if (uniqueNames.count(name) == 0) { //this name hasn't been seen yet
                    uniqueNames.insert(name);
                    wroteSomething = true;
                    
                    out << name << '\t' << group << endl;
                    selectedCount++;
                    
                    if (m->debug) {  sanity["group"].insert(name); }
                }else {
                    m->mothurOut("[WARNING]: " + name + " is in your group file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                }
            }
					
			m->gobble(in);
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) { m->mothurOut("Your file does not contain any sequence from the .accnos file."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName);  outputTypes["group"].push_back(outputFileName);
		
		m->mothurOut("Selected " + toString(selectedCount) + " sequences from your group file."); m->mothurOutEndLine();

		
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "GetSeqsCommand", "readGroup");
		exit(1);
	}
}
//**********************************************************************************************************************
int GetSeqsCommand::readTax(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(taxfile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(taxfile));
        variables["[extension]"] = m->getExtension(taxfile);
		string outputFileName = getOutputFileName("taxonomy", variables);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(taxfile, in);
		string name, tax;
		
		bool wroteSomething = false;
		int selectedCount = 0;
        
        if (m->debug) { set<string> temp; sanity["tax"] = temp; }
		
        set<string> uniqueNames;
		while(!in.eof()){

			if (m->control_pressed) { in.close(); out.close(); m->mothurRemove(outputFileName);  return 0; }

			in >> name;				//read from first column
			in >> tax;			//read from second column
            
            if (!dups) {//adjust name if needed
                map<string, string>::iterator it = uniqueMap.find(name);
                if (it != uniqueMap.end()) { name = it->second; }
            }
			
            if (names.count(name) != 0) {
                if (uniqueNames.count(name) == 0) { //this name hasn't been seen yet
                    uniqueNames.insert(name);

                    wroteSomething = true;
				
                    out << name << '\t' << tax << endl;
                    selectedCount++;
                
                    if (m->debug) { sanity["tax"].insert(name); }
                }else {
                    m->mothurOut("[WARNING]: " + name + " is in your taxonomy file more than once.  Mothur requires sequence names to be unique. I will only add it once.\n");
                }
			}
					
			m->gobble(in);
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) { m->mothurOut("Your file does not contain any sequence from the .accnos file."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName);  outputTypes["taxonomy"].push_back(outputFileName);
		
		m->mothurOut("Selected " + toString(selectedCount) + " sequences from your taxonomy file."); m->mothurOutEndLine();
			
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "GetSeqsCommand", "readTax");
		exit(1);
	}
}
//**********************************************************************************************************************
//alignreport file has a column header line then all other lines contain 16 columns.  we just want the first column since that contains the name
int GetSeqsCommand::readAlign(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(alignfile);  }
        map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(alignfile));
		string outputFileName = getOutputFileName("alignreport", variables);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		

		ifstream in;
		m->openInputFile(alignfile, in);
		string name, junk;
		
		bool wroteSomething = false;
		int selectedCount = 0;
		
		//read column headers
		for (int i = 0; i < 16; i++) {  
			if (!in.eof())	{	in >> junk;	 out << junk << '\t';	}
			else			{	break;			}
		}
		out << endl;
		
        set<string> uniqueNames;
		while(!in.eof()){
		
			if (m->control_pressed) { in.close(); out.close(); m->mothurRemove(outputFileName);  return 0; }


			in >> name;				//read from first column
            
            if (!dups) {//adjust name if needed
                map<string, string>::iterator it = uniqueMap.find(name);
                if (it != uniqueMap.end()) { name = it->second; }
            }
			
			//if this name is in the accnos file
            if (names.count(name) != 0) {
                if (uniqueNames.count(name) == 0) { //this name hasn't been seen yet
                    uniqueNames.insert(name);
                    wroteSomething = true;
                    selectedCount++;
                    
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
		
		if (wroteSomething == false) { m->mothurOut("Your file does not contain any sequence from the .accnos file."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName);  outputTypes["alignreport"].push_back(outputFileName);
		
		m->mothurOut("Selected " + toString(selectedCount) + " sequences from your alignreport file."); m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "GetSeqsCommand", "readAlign");
		exit(1);
	}
}
//**********************************************************************************************************************
//just looking at common mistakes. 
int GetSeqsCommand::runSanityCheck(){
	try {
        string thisOutputDir = outputDir;
        if (outputDir == "") {  thisOutputDir += m->hasPath(fastafile);  }
        string filename = outputDir + "get.seqs.debug.report";
        
        ofstream out;
		m->openOutputFile(filename, out); 

       
        //compare fasta, name, qual and taxonomy if given to make sure they contain the same seqs
        if (fastafile != "") {
            if (namefile != "") { //compare with fasta
                if (sanity["fasta"] != sanity["name"]) { //create mismatch file
                    createMisMatchFile(out, fastafile, namefile, sanity["fasta"], sanity["name"]);
                }
            }
            if (qualfile != "") {
                if (sanity["fasta"] != sanity["qual"]) { //create mismatch file
                    createMisMatchFile(out, fastafile, qualfile, sanity["fasta"], sanity["qual"]);
                }
            }
            if (taxfile != "") {
                if (sanity["fasta"] != sanity["tax"]) { //create mismatch file
                    createMisMatchFile(out, fastafile, taxfile, sanity["fasta"], sanity["tax"]);
                }
            }
        }
        
        //compare dupnames, groups and list if given to make sure they match
        if (namefile != "") {
            if (groupfile != "") {
                if (sanity["dupname"] != sanity["group"]) { //create mismatch file
                    createMisMatchFile(out, namefile, groupfile, sanity["dupname"], sanity["group"]);
                } 
            }
            if (listfile != "") {
                if (sanity["dupname"] != sanity["list"]) { //create mismatch file
                    createMisMatchFile(out, namefile, listfile, sanity["dupname"], sanity["list"]);
                } 
            }
        }else{

            if ((groupfile != "") && (fastafile != "")) {
                if (sanity["fasta"] != sanity["group"]) { //create mismatch file
                    createMisMatchFile(out, fastafile, groupfile, sanity["fasta"], sanity["group"]);
                } 
            }
        }
        
        out.close();
        
        if (m->isBlank(filename)) { m->mothurRemove(filename); }
        else { m->mothurOut("\n[DEBUG]: " + filename + " contains the file mismatches.\n");outputNames.push_back(filename); outputTypes["debug"].push_back(filename); }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "GetSeqsCommand", "runSanityCheck");
		exit(1);
	}
}
//**********************************************************************************************************************
//just looking at common mistakes. 
int GetSeqsCommand::createMisMatchFile(ofstream& out, string filename1, string filename2, set<string> set1, set<string> set2){
	try {
        out << "****************************************" << endl << endl;
        out << "Names unique to " << filename1 << ":\n";
        
        //remove names in set1 that are also in set2
        for (set<string>::iterator it = set1.begin(); it != set1.end();) {
            string name = *it;
            
            if (set2.count(name) == 0)  { out << name << endl;  } //name unique to set1
            else                        { set2.erase(name);     } //you are in both so erase 
            set1.erase(it++);
        }
        
        out << "\nNames unique to " << filename2 << ":\n";
        //output results
        for (set<string>::iterator it = set2.begin(); it != set2.end(); it++) {  out << *it << endl;  }
        
        out << "****************************************" << endl << endl;
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "GetSeqsCommand", "runSanityCheck");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetSeqsCommand::compareAccnos(){
	try {
		
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(accnosfile);  }
        map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + m->getRootName(m->getSimpleName(accnosfile));
		string outputFileName = getOutputFileName("accnosreport", variables);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(accnosfile2, in);
		string name;
		
		set<string> namesAccnos2;
		set<string> namesDups;
		set<string> namesAccnos = names;
		
		map<string, int> nameCount;
		
		if (namefile != "") {
			ifstream inName;
			m->openInputFile(namefile, inName);
			
			
			while(!inName.eof()){
				
				if (m->control_pressed) { inName.close(); return 0; }
				
				string thisname, repnames;
				
				inName >> thisname;		m->gobble(inName);		//read from first column
				inName >> repnames;			//read from second column
				
				int num = m->getNumNames(repnames);
				nameCount[thisname] = num;
				
				m->gobble(inName);
			}
			inName.close();	
		}
		
		while(!in.eof()){
			in >> name;
			
			if (namesAccnos.count(name) == 0){ //name unique to accnos2
				int pos = name.find_last_of('_');
				string tempName = name;
				if (pos != string::npos) {  tempName = tempName.substr(pos+1); cout << tempName << endl; }
				if (namesAccnos.count(tempName) == 0){
					namesAccnos2.insert(name);
				}else { //you are in both so erase
					namesAccnos.erase(name);
					namesDups.insert(name);
				}
			}else { //you are in both so erase
				namesAccnos.erase(name);
				namesDups.insert(name);
			}
			
			m->gobble(in);
		}
		in.close();	
		
		out << "Names in both files : " + toString(namesDups.size()) << endl;
		m->mothurOut("Names in both files : " + toString(namesDups.size())); m->mothurOutEndLine();
		
		for (set<string>::iterator it = namesDups.begin(); it != namesDups.end(); it++) {
			out << (*it);
			if (namefile != "") { out << '\t' << nameCount[(*it)]; }
			out << endl;
		}
		
		out << "Names unique to " + accnosfile + " : " + toString(namesAccnos.size()) << endl;
		m->mothurOut("Names unique to " + accnosfile + " : " + toString(namesAccnos.size())); m->mothurOutEndLine();
		
		for (set<string>::iterator it = namesAccnos.begin(); it != namesAccnos.end(); it++) {
			out << (*it);
			if (namefile != "") { out << '\t' << nameCount[(*it)]; }
			out << endl;
		}
		
		out << "Names unique to " + accnosfile2 + " : " + toString(namesAccnos2.size()) << endl;
		m->mothurOut("Names unique to " + accnosfile2 + " : " + toString(namesAccnos2.size())); m->mothurOutEndLine();
		
		for (set<string>::iterator it = namesAccnos2.begin(); it != namesAccnos2.end(); it++) {
			out << (*it);
			if (namefile != "") { out << '\t' << nameCount[(*it)]; }
			out << endl;
		}

		out.close(); 
		
		outputNames.push_back(outputFileName);  outputTypes["accnosreport"].push_back(outputFileName);
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "GetSeqsCommand", "compareAccnos");
		exit(1);
	}
}


//**********************************************************************************************************************

