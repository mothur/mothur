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

//**********************************************************************************************************************
vector<string> GetSeqsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(pfasta);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "FNGLT", "none",false,false); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "FNGLT", "none",false,false); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "FNGLT", "none",false,false); parameters.push_back(pgroup);
		CommandParameter plist("list", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(plist);
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(ptaxonomy);
		CommandParameter palignreport("alignreport", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(palignreport);
		CommandParameter pqfile("qfile", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(pqfile);
		CommandParameter paccnos("accnos", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(paccnos);
		CommandParameter pdups("dups", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(pdups);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		CommandParameter paccnos2("accnos2", "InputTypes", "", "", "none", "none", "none",false,true); parameters.push_back(paccnos2);

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
		helpString += "The get.seqs command reads an .accnos file and any of the following file types: fasta, name, group, count, list, taxonomy, quality or alignreport file.\n";
		helpString += "It outputs a file containing only the sequences in the .accnos file.\n";
		helpString += "The get.seqs command parameters are accnos, fasta, name, group, list, taxonomy, qfile, alignreport and dups.  You must provide accnos unless you have a valid current accnos file, and at least one of the other parameters.\n";
		helpString += "The dups parameter allows you to add the entire line from a name file if you add any name from the line. default=false. \n";
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
string GetSeqsCommand::getOutputFileNameTag(string type, string inputName=""){	
	try {
        string outputFileName = "";
		map<string, vector<string> >::iterator it;
        
        //is this a type this command creates
        it = outputTypes.find(type);
        if (it == outputTypes.end()) {  m->mothurOut("[ERROR]: this command doesn't create a " + type + " output file.\n"); }
        else {
            if (type == "fasta")            {   outputFileName =  "pick" + m->getExtension(inputName);   }
            else if (type == "taxonomy")    {   outputFileName =  "pick" + m->getExtension(inputName);   }
            else if (type == "name")        {   outputFileName =  "pick" + m->getExtension(inputName);   }
            else if (type == "count")       {   outputFileName =  "pick.count.table";   }
            else if (type == "group")       {   outputFileName =  "pick" + m->getExtension(inputName);   }
            else if (type == "list")        {   outputFileName =  "pick" + m->getExtension(inputName);   }
            else if (type == "qfile")       {   outputFileName =  "pick" + m->getExtension(inputName);   }
            else if (type == "accnosreport"){   outputFileName =  "accnos.report";                       }
            else if (type == "alignreport") {   outputFileName =  "pick.align.report";                   }
            else { m->mothurOut("[ERROR]: No definition for type " + type + " output file tag.\n"); m->control_pressed = true;  }
        }
        return outputFileName;
	}
	catch(exception& e) {
		m->errorOut(e, "GetSeqsCommand", "getOutputFileNameTag");
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
			
			if ((fastafile == "") && (namefile == "") && (groupfile == "") && (alignfile == "") && (listfile == "") && (taxfile == "") && (qualfile == "") && (accnosfile2 == "") && (countfile == ""))  { m->mothurOut("You must provide one of the following: fasta, name, group, count, alignreport, taxonomy, quality or listfile."); m->mothurOutEndLine(); abort = true; }
            
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
int GetSeqsCommand::readFasta(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(fastafile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(fastafile)) + getOutputFileNameTag("fasta", fastafile);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		
		ifstream in;
		m->openInputFile(fastafile, in);
		string name;
		
		bool wroteSomething = false;
		int selectedCount = 0;
        
        if (m->debug) { set<string> temp; sanity["fasta"] = temp; }
		
		while(!in.eof()){
		
			if (m->control_pressed) { in.close(); out.close(); m->mothurRemove(outputFileName);  return 0; }
			
			Sequence currSeq(in);
			name = currSeq.getName();
			
			if (name != "") {
				//if this name is in the accnos file
				if (names.count(name) != 0) {
					wroteSomething = true;
					
					currSeq.printSequence(out);
					selectedCount++;
                    
                    if (m->debug) { sanity["fasta"].insert(name); }
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
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(qualfile)) + getOutputFileNameTag("qfile", qualfile);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		
		ifstream in;
		m->openInputFile(qualfile, in);
		string name;
		
		bool wroteSomething = false;
		int selectedCount = 0;
		
        if (m->debug) { set<string> temp; sanity["qual"] = temp; }
		
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
			
			if (names.count(saveName) != 0) {
				wroteSomething = true;
						
				out << name << endl << scores;
				selectedCount++;
                if (m->debug) { sanity["qual"].insert(name); }
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
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(countfile)) + getOutputFileNameTag("count", countfile);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(countfile, in);
		
		bool wroteSomething = false;
		int selectedCount = 0;
		
        string headers = m->getline(in); m->gobble(in);
        out << headers << endl;
        
        string name, rest; int thisTotal;
        while (!in.eof()) {
            
            if (m->control_pressed) { in.close();  out.close();  m->mothurRemove(outputFileName);  return 0; }
            
            in >> name; m->gobble(in); 
            in >> thisTotal; m->gobble(in);
            rest = m->getline(in); m->gobble(in);
            if (m->debug) { m->mothurOut("[DEBUG]: " + name + '\t' + rest + "\n"); }
            
            if (names.count(name) != 0) {
                out << name << '\t' << thisTotal << '\t' << rest << endl;
                wroteSomething = true;
                selectedCount+= thisTotal;
            }
        }
        in.close();
		out.close();
		
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
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(listfile)) + getOutputFileNameTag("list", listfile);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(listfile, in);
		
		bool wroteSomething = false;
		int selectedCount = 0;
        
        if (m->debug) { set<string> temp; sanity["list"] = temp; }
		
		while(!in.eof()){
			
			selectedCount = 0;
			
			if (m->control_pressed) { in.close(); out.close(); m->mothurRemove(outputFileName);  return 0; }

			//read in list vector
			ListVector list(in);
			
			//make a new list vector
			ListVector newList;
			newList.setLabel(list.getLabel());
			
			//for each bin
			for (int i = 0; i < list.getNumBins(); i++) {
			
				//parse out names that are in accnos file
				string binnames = list.get(i);
				
				string newNames = "";
				while (binnames.find_first_of(',') != -1) { 
					string name = binnames.substr(0,binnames.find_first_of(','));
					binnames = binnames.substr(binnames.find_first_of(',')+1, binnames.length());
					
					//if that name is in the .accnos file, add it
					if (names.count(name) != 0) {  newNames += name + ",";  selectedCount++; if (m->debug) { sanity["list"].insert(name); } }
				}
			
				//get last name
				if (names.count(binnames) != 0) {  newNames += binnames + ",";  selectedCount++;  if (m->debug) { sanity["list"].insert(binnames); } }

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
		
		if (wroteSomething == false) { m->mothurOut("Your file does not contain any sequence from the .accnos file."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName); outputTypes["list"].push_back(outputFileName);
		
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
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(namefile)) + getOutputFileNameTag("name", namefile);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		

		ifstream in;
		m->openInputFile(namefile, in);
		string name, firstCol, secondCol;
		
		bool wroteSomething = false;
		int selectedCount = 0;
        
        if (m->debug) { set<string> temp; sanity["name"] = temp; }
        if (m->debug) { set<string> temp; sanity["dupname"] = temp; }
		
		while(!in.eof()){
		
			if (m->control_pressed) { in.close(); out.close(); m->mothurRemove(outputFileName);  return 0; }

			in >> firstCol;				
			in >> secondCol;
			
			string hold = "";
			if (dups) { hold = secondCol; }
			
			vector<string> parsedNames;
			m->splitAtComma(secondCol, parsedNames);
			
			vector<string> validSecond;
			for (int i = 0; i < parsedNames.size(); i++) {
				if (names.count(parsedNames[i]) != 0) {
					validSecond.push_back(parsedNames[i]);
                    if (m->debug) { sanity["dupname"].insert(parsedNames[i]); }
				}
			}

			if ((dups) && (validSecond.size() != 0)) { //dups = true and we want to add someone, then add everyone
				for (int i = 0; i < parsedNames.size(); i++) {  names.insert(parsedNames[i]); if (m->debug) { sanity["dupname"].insert(parsedNames[i]); } }
				out << firstCol << '\t' << hold << endl;
				wroteSomething = true;
				selectedCount += parsedNames.size();
                if (m->debug) { sanity["name"].insert(firstCol); }
			}else {
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
					
						//you know you have at least one valid second since first column is valid
						for (int i = 0; i < validSecond.size()-1; i++) {  out << validSecond[i] << ',';  }
						out << validSecond[validSecond.size()-1] << endl;
                        
                        if (m->debug) { sanity["name"].insert(validSecond[0]); }
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
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(groupfile)) + getOutputFileNameTag("group", groupfile);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		

		ifstream in;
		m->openInputFile(groupfile, in);
		string name, group;
		
		bool wroteSomething = false;
		int selectedCount = 0;
        
        if (m->debug) { set<string> temp; sanity["group"] = temp; }
		
		while(!in.eof()){

			if (m->control_pressed) { in.close(); out.close(); m->mothurRemove(outputFileName);  return 0; }


			in >> name;				//read from first column
			in >> group;			//read from second column
			
			//if this name is in the accnos file
			if (names.count(name) != 0) {
				wroteSomething = true;
				
				out << name << '\t' << group << endl;
				selectedCount++;
                
                if (m->debug) {  sanity["group"].insert(name); }
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
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(taxfile)) + getOutputFileNameTag("taxonomy", taxfile);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(taxfile, in);
		string name, tax;
		
		bool wroteSomething = false;
		int selectedCount = 0;
        
        if (m->debug) { set<string> temp; sanity["tax"] = temp; }
		
		while(!in.eof()){

			if (m->control_pressed) { in.close(); out.close(); m->mothurRemove(outputFileName);  return 0; }

			in >> name;				//read from first column
			in >> tax;			//read from second column
			
			//if this name is in the accnos file
			if (names.count(name) != 0) {
				wroteSomething = true;
				
				out << name << '\t' << tax << endl;
				selectedCount++;
                
                if (m->debug) { sanity["tax"].insert(name); }
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
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(alignfile)) + getOutputFileNameTag("alignreport");
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
		
		while(!in.eof()){
		
			if (m->control_pressed) { in.close(); out.close(); m->mothurRemove(outputFileName);  return 0; }


			in >> name;				//read from first column
			
			//if this name is in the accnos file
			if (names.count(name) != 0) {
				wroteSomething = true;
				selectedCount++;
				
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
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(accnosfile)) + getOutputFileNameTag("accnosreport");
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

