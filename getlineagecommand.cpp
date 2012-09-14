/*
 *  getlineagecommand.cpp
 *  Mothur
 *
 *  Created by westcott on 9/24/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "getlineagecommand.h"
#include "sequence.hpp"
#include "listvector.hpp"

//**********************************************************************************************************************
vector<string> GetLineageCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(pfasta);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "FNGLT", "none",false,false); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "FNGLT", "none",false,false); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "FNGLT", "none",false,false); parameters.push_back(pgroup);
		CommandParameter plist("list", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(plist);
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "FNGLT", "none",false,true); parameters.push_back(ptaxonomy);
		CommandParameter palignreport("alignreport", "InputTypes", "", "", "none", "FNGLT", "none",false,false); parameters.push_back(palignreport);
		CommandParameter ptaxon("taxon", "String", "", "", "", "", "",false,true); parameters.push_back(ptaxon);
		CommandParameter pdups("dups", "Boolean", "", "T", "", "", "",false,false); parameters.push_back(pdups);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "GetLineageCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetLineageCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The get.lineage command reads a taxonomy file and any of the following file types: fasta, name, group, count, list or alignreport file.\n";
		helpString += "It outputs a file containing only the sequences from the taxonomy file that are from the taxon requested.\n";
		helpString += "The get.lineage command parameters are taxon, fasta, name, group, count, list, taxonomy, alignreport and dups.  You must provide taxonomy unless you have a valid current taxonomy file.\n";
		helpString += "The dups parameter allows you to add the entire line from a name file if you add any name from the line. default=false. \n";
		helpString += "The taxon parameter allows you to select the taxons you would like to get and is required.\n";
		helpString += "You may enter your taxons with confidence scores, doing so will get only those sequences that belong to the taxonomy and whose cofidence scores is above the scores you give.\n";
		helpString += "If they belong to the taxonomy and have confidences below those you provide the sequence will not be selected.\n";
		helpString += "The get.lineage command should be in the following format: get.lineage(taxonomy=yourTaxonomyFile, taxon=yourTaxons).\n";
		helpString += "Example get.lineage(taxonomy=amazon.silva.taxonomy, taxon=Bacteria;Firmicutes;Bacilli;Lactobacillales;).\n";
		helpString += "Note: If you are running mothur in script mode you must wrap the taxon in ' characters so mothur will ignore the ; in the taxon.\n";
		helpString += "Example get.lineage(taxonomy=amazon.silva.taxonomy, taxon='Bacteria;Firmicutes;Bacilli;Lactobacillales;').\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "GetLineageCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetLineageCommand::getOutputFileNameTag(string type, string inputName=""){	
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
            else if (type == "count")       {   outputFileName =  "pick.count.table";                    }  
            else if (type == "group")       {   outputFileName =  "pick" + m->getExtension(inputName);   }
            else if (type == "list")        {   outputFileName =  "pick" + m->getExtension(inputName);   }
            else if (type == "alignreport") {   outputFileName =  "pick.align.report";   }
            else { m->mothurOut("[ERROR]: No definition for type " + type + " output file tag.\n"); m->control_pressed = true;  }
        }
        return outputFileName;
	}
	catch(exception& e) {
		m->errorOut(e, "GetLineageCommand", "getOutputFileNameTag");
		exit(1);
	}
}
//**********************************************************************************************************************
GetLineageCommand::GetLineageCommand(){	
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
        outputTypes["count"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "GetLineageCommand", "GetLineageCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
GetLineageCommand::GetLineageCommand(string option)  {
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
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
			}

			
			//check for required parameters			
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
			else if (taxfile == "not found") {  				
				taxfile = m->getTaxonomyFile(); 
				if (taxfile != "") { m->mothurOut("Using " + taxfile + " as input file for the taxonomy parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current taxonomy file and the taxonomy parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setTaxonomyFile(taxfile); }
			
			string usedDups = "true";
			string temp = validParameter.validFile(parameters, "dups", false);	
			if (temp == "not found") { 
				if (namefile != "") {  temp = "true";					}
				else				{  temp = "false"; usedDups = "";	}
			}
			dups = m->isTrue(temp);
            
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
			
			taxons = validParameter.validFile(parameters, "taxon", false);	
			if (taxons == "not found") { taxons = "";  m->mothurOut("No taxons given, please correct."); m->mothurOutEndLine();  abort = true;  }
			else { 
				//rip off quotes
				if (taxons[0] == '\'') {  taxons = taxons.substr(1); }
				if (taxons[(taxons.length()-1)] == '\'') {  taxons = taxons.substr(0, (taxons.length()-1)); }
			}
			m->splitAtChar(taxons, listOfTaxons, '-');
			
			if ((fastafile == "") && (namefile == "") && (groupfile == "") && (alignfile == "") && (listfile == "") && (taxfile == "") && (countfile == ""))  { m->mothurOut("You must provide one of the following: fasta, name, group, count, alignreport, taxonomy or listfile."); m->mothurOutEndLine(); abort = true; }
		
            if (countfile == "") {
                if ((namefile == "") && ((fastafile != "") || (taxfile != ""))){
                    vector<string> files; files.push_back(fastafile); files.push_back(taxfile);
                    parser.getNameFile(files);
                }
            }
		}

	}
	catch(exception& e) {
		m->errorOut(e, "GetLineageCommand", "GetLineageCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int GetLineageCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		if (m->control_pressed) { return 0; }
        
        if (countfile != "") {
            if ((fastafile != "") || (listfile != "") || (taxfile != "")) { 
                m->mothurOut("\n[NOTE]: The count file should contain only unique names, so mothur assumes your fasta, list and taxonomy files also contain only uniques.\n\n");
            }
        }
		
		//read through the correct file and output lines you want to keep
		if (taxfile != "")			{		readTax();		}  //fills the set of names to get
		if (namefile != "")			{		readName();		}
		if (fastafile != "")		{		readFasta();	}
        if (countfile != "")		{		readCount();	}
		if (groupfile != "")		{		readGroup();	}
		if (alignfile != "")		{		readAlign();	}
		if (listfile != "")			{		readList();		}
		
		
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
			
            itTypes = outputTypes.find("count");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setCountTableFile(current); }
			}
		}
		
		return 0;		
	}

	catch(exception& e) {
		m->errorOut(e, "GetLineageCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
int GetLineageCommand::readFasta(){
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
		
		while(!in.eof()){
		
			if (m->control_pressed) { in.close(); out.close(); m->mothurRemove(outputFileName);  return 0; }
			
			Sequence currSeq(in);
			name = currSeq.getName();
			
			if (name != "") {
				//if this name is in the accnos file
				if (names.count(name) != 0) {
					wroteSomething = true;
					
					currSeq.printSequence(out);
				}
			}
			m->gobble(in);
		}
		in.close();	
		out.close();
		
		if (wroteSomething == false) { m->mothurOut("Your file contains does not contain any sequences from " + taxons + "."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName);  outputTypes["fasta"].push_back(outputFileName);
		
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "GetLineageCommand", "readFasta");
		exit(1);
	}
}
//**********************************************************************************************************************
int GetLineageCommand::readCount(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(countfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(countfile)) + getOutputFileNameTag("count", countfile);
		
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(countfile, in);
		
		bool wroteSomething = false;
		
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
            }
        }
        in.close();
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains does not contain any sequences from " + taxons + "."); m->mothurOutEndLine();  }
		outputTypes["count"].push_back(outputFileName); outputNames.push_back(outputFileName);
		       
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "GetLineageCommand", "readCount");
		exit(1);
	}
}
//**********************************************************************************************************************
int GetLineageCommand::readList(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(listfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(listfile)) + getOutputFileNameTag("list", listfile);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(listfile, in);
		
		bool wroteSomething = false;
		
		while(!in.eof()){
			
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
					if (names.count(name) != 0) {  newNames += name + ",";  }
				}
			
				//get last name
				if (names.count(binnames) != 0) {  newNames += binnames + ",";  }

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
		
		if (wroteSomething == false) { m->mothurOut("Your file contains does not contain any sequences from " + taxons + "."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName); outputTypes["list"].push_back(outputFileName);
		
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "GetLineageCommand", "readList");
		exit(1);
	}
}
//**********************************************************************************************************************
int GetLineageCommand::readName(){
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
				}
			}

			if ((dups) && (validSecond.size() != 0)) { //dups = true and we want to add someone, then add everyone
				for (int i = 0; i < parsedNames.size(); i++) {  names.insert(parsedNames[i]);  }
				out << firstCol << '\t' << hold << endl;
				wroteSomething = true;
			}else {
				//if the name in the first column is in the set then print it and any other names in second column also in set
				if (names.count(firstCol) != 0) {
				
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
		
		if (wroteSomething == false) { m->mothurOut("Your file contains does not contain any sequences from " + taxons + "."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName);  outputTypes["name"].push_back(outputFileName);
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "GetLineageCommand", "readName");
		exit(1);
	}
}

//**********************************************************************************************************************
int GetLineageCommand::readGroup(){
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
		
		while(!in.eof()){

			if (m->control_pressed) { in.close(); out.close(); m->mothurRemove(outputFileName);  return 0; }


			in >> name;				//read from first column
			in >> group;			//read from second column
			
			//if this name is in the accnos file
			if (names.count(name) != 0) {
				wroteSomething = true;
				
				out << name << '\t' << group << endl;
			}
					
			m->gobble(in);
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) { m->mothurOut("Your file contains does not contain any sequences from " + taxons + "."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName);  outputTypes["group"].push_back(outputFileName);
		
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "GetLineageCommand", "readGroup");
		exit(1);
	}
}
//**********************************************************************************************************************
int GetLineageCommand::readTax(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += m->hasPath(taxfile);  }
		string outputFileName = thisOutputDir + m->getRootName(m->getSimpleName(taxfile)) + getOutputFileNameTag("taxonomy", taxfile);
		ofstream out;
		m->openOutputFile(outputFileName, out);
		
		ifstream in;
		m->openInputFile(taxfile, in);
		string name, tax;
		
		//bool wroteSomething = false;
		vector<bool> taxonsHasConfidence; taxonsHasConfidence.resize(listOfTaxons.size(), false);
		vector< vector< map<string, float> > > searchTaxons; searchTaxons.resize(listOfTaxons.size());
		vector<string> noConfidenceTaxons; noConfidenceTaxons.resize(listOfTaxons.size(), "");
		
		for (int i = 0; i < listOfTaxons.size(); i++) {
			noConfidenceTaxons[i] = listOfTaxons[i];
			int hasConPos = listOfTaxons[i].find_first_of('(');
			if (hasConPos != string::npos) {  
				taxonsHasConfidence[i] = true; 
				searchTaxons[i] = getTaxons(listOfTaxons[i]); 
				noConfidenceTaxons[i] = listOfTaxons[i];
				m->removeConfidences(noConfidenceTaxons[i]);
			}
		}
		
		
		while(!in.eof()){

			if (m->control_pressed) { in.close(); out.close(); m->mothurRemove(outputFileName);  return 0; }

			in >> name;				//read from first column
			in >> tax;			//read from second column
			
			for (int j = 0; j < listOfTaxons.size(); j++) {
							
				string newtax = tax;
			
				//if the users file contains confidence scores we want to ignore them when searching for the taxons, unless the taxon has them
				if (!taxonsHasConfidence[j]) {
					int hasConfidences = tax.find_first_of('(');
					if (hasConfidences != string::npos) { 
						newtax = tax;
						m->removeConfidences(newtax);
					}
				
					int pos = newtax.find(noConfidenceTaxons[j]);
				
					if (pos != string::npos) { //this sequence contains the taxon the user wants
						names.insert(name); 
						out << name << '\t' << tax << endl;
						//since you belong to at least one of the taxons we want you are included so no need to search for other
						break;
					}
				}else{//if listOfTaxons[i] has them and you don't them remove taxons
					int hasConfidences = tax.find_first_of('(');
					if (hasConfidences == string::npos) { 
					
						int pos = newtax.find(noConfidenceTaxons[j]);
					
						if (pos != string::npos) { //this sequence contains the taxon the user wants
							names.insert(name);
							out << name << '\t' << tax << endl;
							//since you belong to at least one of the taxons we want you are included so no need to search for other
							break;
						}
					}else { //both have confidences so we want to make sure the users confidences are greater then or equal to the taxons
						//first remove confidences from both and see if the taxonomy exists
					
						string noNewTax = tax;
						int hasConfidences = tax.find_first_of('(');
						if (hasConfidences != string::npos) { 
							noNewTax = tax;
							m->removeConfidences(noNewTax);
						}
					
						int pos = noNewTax.find(noConfidenceTaxons[j]);
					
						if (pos != string::npos) { //if yes, then are the confidences okay
						
							bool good = true;
							vector< map<string, float> > usersTaxon = getTaxons(newtax);
						
							//the usersTaxon is most likely longer than the searchTaxons, and searchTaxon[0] may relate to userTaxon[4]
							//we want to "line them up", so we will find the the index where the searchstring starts
							int index = 0;
							for (int i = 0; i < usersTaxon.size(); i++) {
							
								if (usersTaxon[i].begin()->first == searchTaxons[j][0].begin()->first) { 
									index = i;  
									int spot = 0;
									bool goodspot = true;
									//is this really the start, or are we dealing with a taxon of the same name?
									while ((spot < searchTaxons[j].size()) && ((i+spot) < usersTaxon.size())) {
										if (usersTaxon[i+spot].begin()->first != searchTaxons[j][spot].begin()->first) { goodspot = false; break; }
										else { spot++; }
									}
								
									if (goodspot) { break; }
								}
							}
						
							for (int i = 0; i < searchTaxons[j].size(); i++) {
							
								if ((i+index) < usersTaxon.size()) { //just in case, should never be false
									if (usersTaxon[i+index].begin()->second < searchTaxons[j][i].begin()->second) { //is the users cutoff less than the search taxons
										good = false;
										break;
									}
								}else {
									good = false;
									break;
								}
							}
						
							//passed the test so add you
							if (good) {
								names.insert(name);
								out << name << '\t' << tax << endl;
								break;
							}
						}
					}
				}
			
			}
			
			m->gobble(in);
		}
		in.close();
		out.close();
		
		if (names.size() == 0) { m->mothurOut("Your taxonomy file does not contain any sequences from " + taxons + "."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName); outputTypes["taxonomy"].push_back(outputFileName);
			
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "GetLineageCommand", "readTax");
		exit(1);
	}
}
/**************************************************************************************************/
vector< map<string, float> > GetLineageCommand::getTaxons(string tax) {
	try {
		
		vector< map<string, float> > t;
		string taxon = "";
		int taxLength = tax.length();
		for(int i=0;i<taxLength;i++){
			if(tax[i] == ';'){
		
				int openParen = taxon.find_first_of('(');
				int closeParen = taxon.find_last_of(')');
				
				string newtaxon, confidence;
				if ((openParen != string::npos) && (closeParen != string::npos)) {
					newtaxon = taxon.substr(0, openParen); //rip off confidence
					confidence = taxon.substr((openParen+1), (closeParen-openParen-1));  
				}else{
					newtaxon = taxon;
					confidence = "0";
				} 
				float con = 0;
				convert(confidence, con);
				
				map<string, float> temp;
				temp[newtaxon] = con;
				t.push_back(temp);
				
				taxon = "";
			}
			else{
				taxon += tax[i];
			}
		}
		
		return t;
	}
	catch(exception& e) {
		m->errorOut(e, "GetLineageCommand", "getTaxons");
		exit(1);
	}
}
//**********************************************************************************************************************
//alignreport file has a column header line then all other lines contain 16 columns.  we just want the first column since that contains the name
int GetLineageCommand::readAlign(){
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
		
		if (wroteSomething == false) { m->mothurOut("Your file contains does not contain any sequences from " + taxons + "."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName); outputTypes["alignreport"].push_back(outputFileName);
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "GetLineageCommand", "readAlign");
		exit(1);
	}
}
//**********************************************************************************************************************




