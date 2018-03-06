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
#include "counttable.h"
#include "inputdata.h"

//**********************************************************************************************************************
vector<string> GetLineageCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "FNGLT", "none","fasta",false,false, true); parameters.push_back(pfasta);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "FNGLT", "none","name",false,false, true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "FNGLT", "none","count",false,false, true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "FNGLT", "none","group",false,false, true); parameters.push_back(pgroup);
		CommandParameter plist("list", "InputTypes", "", "", "none", "FNGLT", "none","list",false,false, true); parameters.push_back(plist);
        CommandParameter pshared("shared", "InputTypes", "", "", "none", "FNGLT", "none","shared",false,false, true); parameters.push_back(pshared);
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "tax", "FNGLT", "none","taxonomy",false,false, true); parameters.push_back(ptaxonomy);
        CommandParameter pconstaxonomy("constaxonomy", "InputTypes", "", "", "tax", "FNGLT", "none","constaxonomy",false,false, true); parameters.push_back(pconstaxonomy);
		CommandParameter palignreport("alignreport", "InputTypes", "", "", "none", "FNGLT", "none","alignreport",false,false); parameters.push_back(palignreport);
        CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter ptaxon("taxon", "String", "", "", "", "", "","",false,true, true); parameters.push_back(ptaxon);
		CommandParameter pdups("dups", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pdups);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
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
		helpString += "The get.lineage command reads a taxonomy or constaxonomy file and any of the following file types: fasta, name, group, count, list, shared or alignreport file. The constaxonomy can only be used with a shared or list file.\n";
		helpString += "It outputs a file containing only the sequences from the taxonomy file that are from the taxon requested.\n";
		helpString += "The get.lineage command parameters are taxon, fasta, name, group, count, list, shared, taxonomy, alignreport, label and dups.  You must provide taxonomy or constaxonomy unless you have a valid current taxonomy file.\n";
		helpString += "The dups parameter allows you to add the entire line from a name file if you add any name from the line. default=false. \n";
		helpString += "The taxon parameter allows you to select the taxons you would like to get and is required.\n";
		helpString += "You may enter your taxons with confidence scores, doing so will get only those sequences that belong to the taxonomy and whose cofidence scores is above the scores you give.\n";
		helpString += "If they belong to the taxonomy and have confidences below those you provide the sequence will not be selected.\n";
         helpString += "The label parameter is used to analyze specific labels in your input. \n";
		helpString += "The get.lineage command should be in the following format: get.lineage(taxonomy=yourTaxonomyFile, taxon=yourTaxons).\n";
		helpString += "Example get.lineage(taxonomy=amazon.silva.taxonomy, taxon=Bacteria;Firmicutes;Bacilli;Lactobacillales;).\n";
		helpString += "Note: If you are running mothur in script mode you must wrap the taxon in ' characters so mothur will ignore the ; in the taxon.\n";
		helpString += "Example get.lineage(taxonomy=amazon.silva.taxonomy, taxon='Bacteria;Firmicutes;Bacilli;Lactobacillales;').\n";
		;
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "GetLineageCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string GetLineageCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta")                {   pattern = "[filename],pick,[extension]";    }
        else if (type == "taxonomy")        {   pattern = "[filename],pick,[extension]";    }
        else if (type == "constaxonomy")    {   pattern = "[filename],pick,[extension]";    }
        else if (type == "name")            {   pattern = "[filename],pick,[extension]";    }
        else if (type == "group")           {   pattern = "[filename],pick,[extension]";    }
        else if (type == "count")           {   pattern = "[filename],pick,[extension]";    }
        else if (type == "list")            {   pattern = "[filename],[distance],pick,[extension]";    }
        else if (type == "shared")          {   pattern = "[filename],[distance],pick,[extension]";    }
        else if (type == "alignreport")     {   pattern = "[filename],pick.align.report";    }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "GetLineageCommand", "getOutputPattern");
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
        outputTypes["constaxonomy"] = tempOutNames;
        outputTypes["shared"] = tempOutNames;
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
            outputTypes["constaxonomy"] = tempOutNames;
            outputTypes["shared"] = tempOutNames;

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
                
                it = parameters.find("constaxonomy");
				//user has given a template file
				if(it != parameters.end()){
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["constaxonomy"] = inputDir + it->second;		}
				}
                
                it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
			}

			
			//check for required parameters			
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
			if (taxfile == "not open") { taxfile = ""; abort = true; }
			else if (taxfile == "not found") {  taxfile = "";		}
			else { current->setTaxonomyFile(taxfile); }
            
            sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }
			else if (sharedfile == "not found") {  sharedfile = "";		}
			else { current->setSharedFile(sharedfile); }

            
            constaxonomy = validParameter.validFile(parameters, "constaxonomy");
			if (constaxonomy == "not open") { constaxonomy = ""; abort = true; }
			else if (constaxonomy == "not found") {  constaxonomy = "";		}
    
            if ((constaxonomy == "") && (taxfile == "")) {
                taxfile = current->getTaxonomyFile();
                if (taxfile != "") { m->mothurOut("Using " + taxfile + " as input file for the taxonomy parameter."); m->mothurOutEndLine(); }
                else {
                    m->mothurOut("You have no current taxonomy file and did not provide a constaxonomy file. The taxonomy or constaxonomy parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
            
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
			
			taxons = validParameter.valid(parameters, "taxon");
			if (taxons == "not found") { taxons = "";  m->mothurOut("No taxons given, please correct."); m->mothurOutEndLine();  abort = true;  }
			else { 
				//rip off quotes
				if (taxons[0] == '\'') {  taxons = taxons.substr(1); }
				if (taxons[(taxons.length()-1)] == '\'') {  taxons = taxons.substr(0, (taxons.length()-1)); }
			}
			util.splitAtChar(taxons, listOfTaxons, '-');
			
			if ((fastafile == "") && (constaxonomy == "") && (namefile == "") && (groupfile == "") && (alignfile == "") && (listfile == "") && (taxfile == "") && (countfile == ""))  { m->mothurOut("You must provide one of the following: fasta, name, group, count, alignreport, taxonomy, constaxonomy, shared or listfile."); m->mothurOutEndLine(); abort = true; }
            
            if ((constaxonomy != "") && ((fastafile != "") || (namefile != "") || (groupfile != "") || (alignfile != "") || (taxfile != "") || (countfile != ""))) {
                m->mothurOut("[ERROR]: can only use constaxonomy file with a list or shared file, aborting.\n");  abort = true;
            }
            
            if ((constaxonomy != "") && (taxfile != "")) {
                m->mothurOut("[ERROR]: Choose only one: taxonomy or constaxonomy, aborting.\n"); abort = true;
            }
            
            if ((sharedfile != "") && (taxfile != "")) {
                m->mothurOut("[ERROR]: sharedfile can only be used with constaxonomy file, aborting.\n"); abort = true;
            }
            
            if ((sharedfile != "") || (listfile != "")) {
                label = validParameter.valid(parameters, "label");
                if (label == "not found") { label = ""; m->mothurOut("[WARNING]: You did not provide a label, I will use the first label in your inputfile."); m->mothurOutEndLine(); }
            }
            
            if (countfile == "") {
                if ((namefile == "") && ((fastafile != "") || (taxfile != ""))){
                    vector<string> files; files.push_back(fastafile); files.push_back(taxfile);
                    if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
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
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		if (m->getControl_pressed()) { return 0; }
        
        if (countfile != "") {
            if ((fastafile != "") || (listfile != "") || (taxfile != "")) { 
                m->mothurOut("\n[NOTE]: The count file should contain only unique names, so mothur assumes your fasta, list and taxonomy files also contain only uniques.\n\n");
            }
        }
        
		//read through the correct file and output lines you want to keep
		if (taxfile != "")			{
            readTax(); //fills the set of names to get
            if (namefile != "")			{		readName();		}
            if (fastafile != "")		{		readFasta();	}
            if (countfile != "")		{		readCount();	}
            if (groupfile != "")		{		readGroup();	}
            if (alignfile != "")		{		readAlign();	}
            if (listfile != "")			{		readList();		}

        }else {
            readConsTax();
            if (listfile != "")			{		readConsList();		}
            if (sharedfile != "")		{		readShared();		}
        }
				
		
		if (m->getControl_pressed()) { outputTypes.clear(); for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);  } return 0; }
		
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
            
            itTypes = outputTypes.find("shared");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setSharedFile(currentName); }
			}
			
			itTypes = outputTypes.find("taxonomy");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setTaxonomyFile(currentName); }
			}
			
            itTypes = outputTypes.find("count");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setCountFile(currentName); }
			}
            
            //set constaxonomy file as new current constaxonomyfile
            itTypes = outputTypes.find("constaxonomy");
            if (itTypes != outputTypes.end()) {
                if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setConsTaxonomyFile(currentName); }
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
		
		while(!in.eof()){
		
			if (m->getControl_pressed()) { in.close(); out.close(); util.mothurRemove(outputFileName);  return 0; }
			
			Sequence currSeq(in);
			name = currSeq.getName();
			
			if (name != "") {
				//if this name is in the accnos file
				if (names.count(name) != 0) {
					wroteSomething = true;
					
					currSeq.printSequence(out);
				}
			}
			util.gobble(in);
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
		if (outputDir == "") {  thisOutputDir += util.hasPath(countfile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(countfile));
        variables["[extension]"] = util.getExtension(countfile);
		string outputFileName = getOutputFileName("count", variables);
		
		ofstream out;
		util.openOutputFile(outputFileName, out);
		
		ifstream in;
		util.openInputFile(countfile, in);
		
		bool wroteSomething = false;
		
        string headers = util.getline(in); util.gobble(in);
        out << headers << endl;
        string test = headers; vector<string> pieces = util.splitWhiteSpace(test);
        
        string name, rest; int thisTotal; rest = "";
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return 0; }
            
            in >> name; util.gobble(in); 
            in >> thisTotal; util.gobble(in);
            if (pieces.size() > 2) {  rest = util.getline(in); util.gobble(in);  }
            if (m->getDebug()) { m->mothurOut("[DEBUG]: " + name + '\t' + rest + "\n"); }
            
            if (names.count(name) != 0) {
                out << name << '\t' << thisTotal << '\t' << rest << endl;
                wroteSomething = true;
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
		if (outputDir == "") {  thisOutputDir += util.hasPath(listfile);  }
		map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(listfile));
        variables["[extension]"] = util.getExtension(listfile);
        InputData input(listfile, "list", nullVector);
        ListVector* list = input.getListVector();
        
		bool wroteSomething = false;
		
        while(list != NULL) {
			
			//make a new list vector
			ListVector newList;
			newList.setLabel(list->getLabel());
            
            variables["[distance]"] = list->getLabel();
            string outputFileName = getOutputFileName("list", variables);
			
			ofstream out;
			util.openOutputFile(outputFileName, out);
			outputTypes["list"].push_back(outputFileName);  outputNames.push_back(outputFileName);
            
            if (m->getControl_pressed()) {  out.close(); return 0; }
            
            vector<string> binLabels = list->getLabels();
            vector<string> newBinLabels;
			
			//for each bin
			for (int i = 0; i < list->getNumBins(); i++) {
			
				//parse out names that are in accnos file
				string binnames = list->get(i);
                vector<string> bnames;
                util.splitAtComma(binnames, bnames);
				
				string newNames = "";
                for (int j = 0; j < bnames.size(); j++) {
					string name = bnames[j];
					//if that name is in the .accnos file, add it
					if (names.count(name) != 0) {  newNames += name + ",";  }
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
		
		
		if (wroteSomething == false) { m->mothurOut("Your file contains does not contain any sequences from " + taxons + "."); m->mothurOutEndLine();  }
		
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "GetLineageCommand", "readList");
		exit(1);
	}
}
//**********************************************************************************************************************
int GetLineageCommand::readConsList(){
	try {
		getListVector();
        
        if (m->getControl_pressed()) { delete list; return 0;}
        
        ListVector newList;
        newList.setLabel(list->getLabel());
        int selectedCount = 0;
        bool wroteSomething = false;
        string snumBins = toString(list->getNumBins());
        
        vector<string> binLabels = list->getLabels();
        vector<string> newBinLabels;
        for (int i = 0; i < list->getNumBins(); i++) {
            
            if (m->getControl_pressed()) { delete list; return 0;}
            
            //create a label for this otu
            string otuLabel = "Otu";
            string sbinNumber = toString(i+1);
            if (sbinNumber.length() < snumBins.length()) {
                int diff = snumBins.length() - sbinNumber.length();
                for (int h = 0; h < diff; h++) { otuLabel += "0"; }
            }
            otuLabel += sbinNumber;
            
            if (names.count(util.getSimpleLabel(otuLabel)) != 0) {
				selectedCount++;
                newList.push_back(list->get(i));
                newBinLabels.push_back(binLabels[i]);
            }
        }
        
        string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(listfile);  }
        map<string, string> variables;
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(listfile));
        variables["[extension]"] = util.getExtension(listfile);
        variables["[distance]"] = list->getLabel();
		string outputFileName = getOutputFileName("list", variables);
		ofstream out;
		util.openOutputFile(outputFileName, out);
        
		delete list;
        //print new listvector
        if (newList.getNumBins() != 0) {
            wroteSomething = true;
            newList.setLabels(newBinLabels);
            newList.print(out);
        }
		out.close();
		
		if (wroteSomething == false) { m->mothurOut("Your file does not contain OTUs from " + taxons + "."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName); outputTypes["list"].push_back(outputFileName);
		
		m->mothurOut("Selected " + toString(selectedCount) + " OTUs from your list file."); m->mothurOutEndLine();
        
		return 0;
        
	}
	catch(exception& e) {
		m->errorOut(e, "GetLineageCommand", "readConsList");
		exit(1);
	}
}
//**********************************************************************************************************************
int GetLineageCommand::getListVector(){
	try {
		InputData input(listfile, "list", nullVector);
		list = input.getListVector();
		string lastLabel = list->getLabel();
		
		if (label == "") { label = lastLabel;  return 0; }
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> labels; labels.insert(label);
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((list != NULL) && (userLabels.size() != 0)) {
			if (m->getControl_pressed()) {  return 0;  }
			
			if(labels.count(list->getLabel()) == 1){
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
				break;
			}
			
			if ((util.anyLabelsToProcess(list->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = list->getLabel();
				
				delete list;
				list = input.getListVector(lastLabel);
				
				processedLabels.insert(list->getLabel());
				userLabels.erase(list->getLabel());
				
				//restore real lastlabel to save below
				list->setLabel(saveLabel);
				break;
			}
			
			lastLabel = list->getLabel();
			
			//get next line to process
			//prevent memory leak
			delete list;
			list = input.getListVector();
		}
		
		
		if (m->getControl_pressed()) {  return 0;  }
		
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
		
		//run last label if you need to
		if (needToRun )  {
			delete list;
			list = input.getListVector(lastLabel);
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "GetLineageCommand", "getListVector");
		exit(1);
	}
}

//**********************************************************************************************************************
int GetLineageCommand::readShared(){
	try {
        
        SharedRAbundVectors* lookup = getShared();
        
        if (m->getControl_pressed()) { delete lookup; }
        
        vector<string> newLabels;
        
        bool wroteSomething = false;
        int numSelected = 0;
        
        for (int i = 0; i < lookup->getNumBins();) {
            
            if (m->getControl_pressed()) {  delete lookup; return 0; }
            
            //is this otu on the list
            if (names.count(util.getSimpleLabel(lookup->getOTUName(i))) != 0) {
                numSelected++; wroteSomething = true; ++i;
            }else { lookup->removeOTU(i);  }
        }
        
        string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(sharedfile);  }
        map<string, string> variables;
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(sharedfile));
        variables["[extension]"] = util.getExtension(sharedfile);
        variables["[distance]"] = lookup->getLabel();
		string outputFileName = getOutputFileName("shared", variables);
        ofstream out;
		util.openOutputFile(outputFileName, out);
		outputTypes["shared"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        
        lookup->print(out);
		out.close();
        
        delete lookup;
        
        if (wroteSomething == false) { m->mothurOut("Your file does not contain OTUs from " + taxons + "."); m->mothurOutEndLine();  }
        
		m->mothurOut("Selected " + toString(numSelected) + " OTUs from your shared file."); m->mothurOutEndLine();
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "GetLineageCommand", "readShared");
		exit(1);
	}
}
//**********************************************************************************************************************
SharedRAbundVectors* GetLineageCommand::getShared(){
	try {
		InputData input(sharedfile, "sharedfile", nullVector);
		SharedRAbundVectors* lookup = input.getSharedRAbundVectors();
		string lastLabel = lookup->getLabel();
		
		if (label == "") { label = lastLabel;  return lookup; }
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> labels; labels.insert(label);
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((lookup != NULL) && (userLabels.size() != 0)) {
			if (m->getControl_pressed()) {   return 0;  }
			
			if(labels.count(lookup->getLabel()) == 1){
				processedLabels.insert(lookup->getLabel()); userLabels.erase(lookup->getLabel());
				break;
			}
			
			if ((util.anyLabelsToProcess(lookup->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = lookup->getLabel();
				
				delete lookup;
				lookup = input.getSharedRAbundVectors(lastLabel);
				
				processedLabels.insert(lookup->getLabel()); userLabels.erase(lookup->getLabel());
				
				//restore real lastlabel to save below
				lookup->setLabels(saveLabel);
				break;
			}
			
			lastLabel = lookup->getLabel();
			
			//get next line to process
			//prevent memory leak
            delete lookup;
			lookup = input.getSharedRAbundVectors();
		}
		
		
		if (m->getControl_pressed()) {  return 0;  }
		
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
		
		//run last label if you need to
		if (needToRun )  {
            delete lookup;
			lookup = input.getSharedRAbundVectors(lastLabel);
		}
		
		return lookup;
	}
	catch(exception& e) {
		m->errorOut(e, "GetLineageCommand", "getShared");
		exit(1);
	}
}

//**********************************************************************************************************************
int GetLineageCommand::readName(){
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
		
		
		while(!in.eof()){
		
			if (m->getControl_pressed()) { in.close(); out.close(); util.mothurRemove(outputFileName);  return 0; }

			in >> firstCol;				
			in >> secondCol;
			
			string hold = "";
			if (dups) { hold = secondCol; }
			
			vector<string> parsedNames;
			util.splitAtComma(secondCol, parsedNames);
			
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
			util.gobble(in);
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
		
		while(!in.eof()){

			if (m->getControl_pressed()) { in.close(); out.close(); util.mothurRemove(outputFileName);  return 0; }


			in >> name;				//read from first column
			in >> group;			//read from second column
			
			//if this name is in the accnos file
			if (names.count(name) != 0) {
				wroteSomething = true;
				
				out << name << '\t' << group << endl;
			}
					
			util.gobble(in);
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
		
		//bool wroteSomething = false;
		vector<bool> taxonsHasConfidence; taxonsHasConfidence.resize(listOfTaxons.size(), false);
		vector< vector< map<string, float> > > searchTaxons; searchTaxons.resize(listOfTaxons.size());
		vector<string> noConfidenceTaxons; noConfidenceTaxons.resize(listOfTaxons.size(), "");
		
		for (int i = 0; i < listOfTaxons.size(); i++) {
			noConfidenceTaxons[i] = listOfTaxons[i];
            bool hasCon = false;
            searchTaxons[i] = util.getTaxons(listOfTaxons[i], hasCon);  taxonsHasConfidence[i] = hasCon;
            noConfidenceTaxons[i] = listOfTaxons[i];
            if (hasCon) { util.removeConfidences(noConfidenceTaxons[i]); }
		}
		
		while(!in.eof()){

			if (m->getControl_pressed()) { in.close(); out.close(); util.mothurRemove(outputFileName);  return 0; }

            in >> name; util.gobble(in);
            tax = util.getline(in); util.gobble(in);
			
            string noQuotesTax = util.removeQuotes(tax);
            
			if (util.searchTax(noQuotesTax, listOfTaxons, taxonsHasConfidence, noConfidenceTaxons, searchTaxons)) {
                names.insert(name); out << name << '\t' << tax << endl;
            }
        }
		in.close();
		out.close();
		names.insert(name); out << name << '\t' << tax << endl;
		if (names.size() == 0) { m->mothurOut("Your taxonomy file does not contain any sequences from " + taxons + ".\n");  }
		outputNames.push_back(outputFileName); outputTypes["taxonomy"].push_back(outputFileName);
			
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "GetLineageCommand", "readTax");
		exit(1);
	}
}
//**********************************************************************************************************************
int GetLineageCommand::readConsTax(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(constaxonomy);  }
		map<string, string> variables;
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(constaxonomy));
        variables["[extension]"] = util.getExtension(constaxonomy);
		string outputFileName = getOutputFileName("constaxonomy", variables);
		ofstream out;
		util.openOutputFile(outputFileName, out);
		
		ifstream in;
		util.openInputFile(constaxonomy, in);
		string otuLabel, tax;
        int numReps;
        
        //read headers
        string headers = util.getline(in);
        out << headers << endl;
		
		//bool wroteSomething = false;
		vector<bool> taxonsHasConfidence; taxonsHasConfidence.resize(listOfTaxons.size(), false);
		vector< vector< map<string, float> > > searchTaxons; searchTaxons.resize(listOfTaxons.size());
		vector<string> noConfidenceTaxons; noConfidenceTaxons.resize(listOfTaxons.size(), "");
		
		for (int i = 0; i < listOfTaxons.size(); i++) {
            noConfidenceTaxons[i] = listOfTaxons[i];
            bool hasCon = false;
            searchTaxons[i] = util.getTaxons(listOfTaxons[i], hasCon);  taxonsHasConfidence[i] = hasCon;
            noConfidenceTaxons[i] = listOfTaxons[i];
            if (hasCon) { util.removeConfidences(noConfidenceTaxons[i]); }
        }
        
		while(!in.eof()){
            
			if (m->getControl_pressed()) { in.close(); out.close(); util.mothurRemove(outputFileName);  return 0; }
            
			in >> otuLabel;	 		util.gobble(in);
            in >> numReps;          util.gobble(in);
            tax = util.getline(in);   util.gobble(in);
			
            string noQuotesTax = util.removeQuotes(tax);
            
            if (util.searchTax(noQuotesTax, listOfTaxons, taxonsHasConfidence, noConfidenceTaxons, searchTaxons)) {
                names.insert(util.getSimpleLabel(otuLabel));
                out << otuLabel << '\t' << numReps << '\t' << tax << endl;
            }
        }
		in.close();
		out.close();
		
		if (names.size() == 0) { m->mothurOut("Your taxonomy file does not contain any OTUs from " + taxons + "."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName); outputTypes["constaxonomy"].push_back(outputFileName);
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "GetLineageCommand", "readConsTax");
		exit(1);
	}
}
//**********************************************************************************************************************
//alignreport file has a column header line then all other lines contain 16 columns.  we just want the first column since that contains the name
int GetLineageCommand::readAlign(){
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(alignfile);  }
        map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(alignfile));
        variables["[extension]"] = util.getExtension(alignfile);
		string outputFileName = getOutputFileName("alignreport", variables);
		
		ofstream out;
		util.openOutputFile(outputFileName, out);
		

		ifstream in;
		util.openInputFile(alignfile, in);
		string name, junk;
		
		bool wroteSomething = false;
		
		//read column headers
		for (int i = 0; i < 16; i++) {  
			if (!in.eof())	{	in >> junk;	 out << junk << '\t';	}
			else			{	break;			}
		}
		out << endl;
		
		while(!in.eof()){
		
			if (m->getControl_pressed()) { in.close(); out.close(); util.mothurRemove(outputFileName);  return 0; }


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
			
			util.gobble(in);
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




