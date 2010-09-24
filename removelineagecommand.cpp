/*
 *  removelineagecommand.cpp
 *  Mothur
 *
 *  Created by westcott on 9/24/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "removelineagecommand.h"
#include "sequence.hpp"
#include "listvector.hpp"

//**********************************************************************************************************************

RemoveLineageCommand::RemoveLineageCommand(string option)  {
	try {
		abort = false;
				
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta","name", "group", "alignreport", "taxon", "dups", "list","taxonomy","outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
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
			else if (taxfile == "not found") {  taxfile = ""; m->mothurOut("The taxonomy parameter is required for the get.lineage command."); m->mothurOutEndLine();  abort = true; }
			
			string usedDups = "true";
			string temp = validParameter.validFile(parameters, "dups", false);	if (temp == "not found") { temp = "false"; usedDups = ""; }
			dups = m->isTrue(temp);
			
			taxons = validParameter.validFile(parameters, "taxon", false);	
			if (taxons == "not found") { taxons = "";  m->mothurOut("No taxons given, please correct."); m->mothurOutEndLine();  abort = true;  }
			else { 
				//rip off quotes
				if (taxons[0] == '\'') {  taxons = taxons.substr(1); }
				if (taxons[(taxons.length()-1)] == '\'') {  taxons = taxons.substr(0, (taxons.length()-1)); }
			}
			
			
			if ((fastafile == "") && (namefile == "") && (groupfile == "") && (alignfile == "") && (listfile == "") && (taxfile == ""))  { m->mothurOut("You must provide one of the following: fasta, name, group, alignreport, taxonomy or listfile."); m->mothurOutEndLine(); abort = true; }
		
			if ((usedDups != "") && (namefile == "")) {  m->mothurOut("You may only use dups with the name option."); m->mothurOutEndLine();  abort = true; }			

		}

	}
	catch(exception& e) {
		m->errorOut(e, "RemoveLineageCommand", "RemoveLineageCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void RemoveLineageCommand::help(){
	try {
		m->mothurOut("The remove.lineage command reads a taxonomy file and any of the following file types: fasta, name, group, list or alignreport file.\n");
		m->mothurOut("It outputs a file containing only the sequences from the taxonomy file that are not from the taxon you requested to be removed.\n");
		m->mothurOut("The remove.lineage command parameters are taxon, fasta, name, group, list, taxonomy, alignreport and dups.  You must provide taxonomy and taxon.\n");
		m->mothurOut("The dups parameter allows you to add the entire line from a name file if you add any name from the line. default=false. \n");
		m->mothurOut("The taxon parameter allows you to select the taxons you would like to remove.\n");
		m->mothurOut("The remove.lineage command should be in the following format: remove.lineage(taxonomy=yourTaxonomyFile, taxon=yourTaxons).\n");
		m->mothurOut("Example remove.lineage(taxonomy=amazon.silva.taxonomy, taxon=Bacteria;Firmicutes;Bacilli;Lactobacillales;).\n");
		m->mothurOut("Note: If you are running mothur in script mode you must wrap the taxon in ' characters so mothur will ignore the ; in the taxon.\n");
		m->mothurOut("Example remove.lineage(taxonomy=amazon.silva.taxonomy, taxon='Bacteria;Firmicutes;Bacilli;Lactobacillales;').\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveLineageCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

int RemoveLineageCommand::execute(){
	try {
		
		if (abort == true) { return 0; }
		
		if (m->control_pressed) { return 0; }
		
		//read through the correct file and output lines you want to keep
		if (taxfile != "")			{		readTax();		}  //fills the set of names to remove
		if (namefile != "")			{		readName();		}
		if (fastafile != "")		{		readFasta();	}
		if (groupfile != "")		{		readGroup();	}
		if (alignfile != "")		{		readAlign();	}
		if (listfile != "")			{		readList();		}
		
		
		if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	remove(outputNames[i].c_str());  } return 0; }
		
		if (outputNames.size() != 0) {
			m->mothurOutEndLine();
			m->mothurOut("Output File Names: "); m->mothurOutEndLine();
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
			m->mothurOutEndLine();
		}
		
		return 0;		
	}

	catch(exception& e) {
		m->errorOut(e, "RemoveLineageCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
int RemoveLineageCommand::readFasta(){
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
		
		if (wroteSomething == false) {  m->mothurOut("Your fasta file contains only sequences from " + taxons + "."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName); 
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveLineageCommand", "readFasta");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveLineageCommand::readList(){
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
		
		if (wroteSomething == false) {  m->mothurOut("Your list file contains only sequences from " + taxons + "."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName); 
				
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "RemoveLineageCommand", "readList");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveLineageCommand::readName(){
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

			in >> firstCol;				
			in >> secondCol;			

			vector<string> parsedNames;
			//parse second column saving each name
			while (secondCol.find_first_of(',') != -1) { 
				name = secondCol.substr(0,secondCol.find_first_of(','));
				secondCol = secondCol.substr(secondCol.find_first_of(',')+1, secondCol.length());
				parsedNames.push_back(name);
			}
			
			//get name after last ,
			parsedNames.push_back(secondCol);

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

		if (wroteSomething == false) {  m->mothurOut("Your name file contains only sequences from " + taxons + "."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName); 
				
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveLineageCommand", "readName");
		exit(1);
	}
}

//**********************************************************************************************************************
int RemoveLineageCommand::readGroup(){
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
		
		if (wroteSomething == false) {  m->mothurOut("Your group file contains only sequences from " + taxons + "."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName); 
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveLineageCommand", "readGroup");
		exit(1);
	}
}
//**********************************************************************************************************************
int RemoveLineageCommand::readTax(){
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

			if (m->control_pressed) { in.close(); out.close(); remove(outputFileName.c_str());  return 0; }

			in >> name;				//read from first column
			in >> tax;			//read from second column
			
			string newtax = tax;
			
			//if the users file contains confidence scores we want to ignore them when searching for the taxons
			int hasConfidences = tax.find_first_of('(');
			if (hasConfidences != string::npos) { 
				newtax = removeConfidences(tax);
			}
			
			int pos = newtax.find(taxons);
			
			if (pos == string::npos) { 
				wroteSomething = true;
				out << name << '\t' << tax << endl;
			}else{ //this sequence contains the taxon the user wants to remove
				names.insert(name);
			}
			
			m->gobble(in);
		}
		in.close();
		out.close();
		
		if (!wroteSomething) { m->mothurOut("Your taxonomy file contains only sequences from " + taxons + "."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName); 
			
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "RemoveLineageCommand", "readTax");
		exit(1);
	}
}
/**************************************************************************************************/
string RemoveLineageCommand::removeConfidences(string tax) {
	try {
		
		string taxon = "";
		int taxLength = tax.length();
		for(int i=0;i<taxLength;i++){
			if(tax[i] == ';'){
				taxon = taxon.substr(0, taxon.find_first_of('(')); //rip off confidence
				taxon += ";";
			}
			else{
				taxon += tax[i];
			}
		}
				
		return taxon;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveLineageCommand", "removeConfidences");
		exit(1);
	}
}
//**********************************************************************************************************************
//alignreport file has a column header line then all other lines contain 16 columns.  we just want the first column since that contains the name
int RemoveLineageCommand::readAlign(){
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
		
		if (wroteSomething == false) {  m->mothurOut("Your align file contains only sequences from " + taxons + "."); m->mothurOutEndLine();  }
		outputNames.push_back(outputFileName); 
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveLineageCommand", "readAlign");
		exit(1);
	}
}
//**********************************************************************************************************************

