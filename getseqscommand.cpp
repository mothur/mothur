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

GetSeqsCommand::GetSeqsCommand(string option){
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta","name", "group", "alignreport", "accnos", "list"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//check for required parameters
			accnosfile = validParameter.validFile(parameters, "accnos", true);
			if (accnosfile == "not open") { abort = true; }
			else if (accnosfile == "not found") {  accnosfile = "";  mothurOut("You must provide an accnos file."); mothurOutEndLine(); abort = true; }	
			
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
			
			if ((fastafile == "") && (namefile == "") && (groupfile == "") && (alignfile == "") && (listfile == ""))  { mothurOut("You must provide one of the following: fasta, name, group, alignreport or listfile."); mothurOutEndLine(); abort = true; }
			
			if (parameters.size() > 2) { mothurOut("You may only enter one of the following: fasta, name, group, alignreport or listfile."); mothurOutEndLine(); abort = true;  }
		}

	}
	catch(exception& e) {
		errorOut(e, "GetSeqsCommand", "GetSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void GetSeqsCommand::help(){
	try {
		mothurOut("The get.seqs command reads an .accnos file and one of the following file types: fasta, name, group, list or alignreport file.\n");
		mothurOut("It outputs a file containing only the sequences in the .accnos file.\n");
		mothurOut("The get.seqs command parameters are accnos, fasta, name, group, list and alignreport.  You must provide accnos and one of the other parameters.\n");
		mothurOut("The get.seqs command should be in the following format: get.seqs(accnos=yourAccnos, fasta=yourFasta).\n");
		mothurOut("Example get.seqs(accnos=amazon.accnos, fasta=amazon.fasta).\n");
		mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n\n");
	}
	catch(exception& e) {
		errorOut(e, "GetSeqsCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

int GetSeqsCommand::execute(){
	try {
		
		if (abort == true) { return 0; }
		
		//get names you want to keep
		readAccnos();
		
		//read through the correct file and output lines you want to keep
		if (fastafile != "")		{		readFasta();	}
		else if (namefile != "")	{		readName();		}
		else if (groupfile != "")	{		readGroup();	}
		else if (alignfile != "")	{		readAlign();	}
		else if (listfile != "")	{		readList();		}
		
		return 0;		
	}

	catch(exception& e) {
		errorOut(e, "GetSeqsCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
void GetSeqsCommand::readFasta(){
	try {
		string outputFileName = getRootName(fastafile) + "pick" +  getExtension(fastafile);
		ofstream out;
		openOutputFile(outputFileName, out);
		
		ifstream in;
		openInputFile(fastafile, in);
		string name;
		
		bool wroteSomething = false;
		
		while(!in.eof()){
			Sequence currSeq(in);
			name = currSeq.getName();
			
			if (name != "") {
				//if this name is in the accnos file
				if (names.count(name) == 1) {
					wroteSomething = true;
					
					currSeq.printSequence(out);
					
					names.erase(name);
				}
			}
			gobble(in);
		}
		in.close();	
		out.close();
		
		if (wroteSomething == false) {
			mothurOut("Your file does not contain any sequence from the .accnos file."); mothurOutEndLine();
			remove(outputFileName.c_str()); 
		}

	}
	catch(exception& e) {
		errorOut(e, "GetSeqsCommand", "readFasta");
		exit(1);
	}
}
//**********************************************************************************************************************
void GetSeqsCommand::readList(){
	try {
		string outputFileName = getRootName(listfile) + "pick" +  getExtension(listfile);
		ofstream out;
		openOutputFile(outputFileName, out);
		
		ifstream in;
		openInputFile(listfile, in);
		
		bool wroteSomething = false;
		
		while(!in.eof()){
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
					if (names.count(name) == 1) {  newNames += name + ",";  }
				}
			
				//get last name
				if (names.count(binnames) == 1) {  newNames += binnames;  }

				//if there are names in this bin add to new list
				if (newNames != "") {  newList.push_back(newNames);	}
			}
				
			//print new listvector
			if (newList.getNumBins() != 0) {
				wroteSomething = true;
				newList.print(out);
			}
			
			gobble(in);
		}
		in.close();	
		out.close();
		
		if (wroteSomething == false) {
			mothurOut("Your file does not contain any sequence from the .accnos file."); mothurOutEndLine();
			remove(outputFileName.c_str()); 
		}

	}
	catch(exception& e) {
		errorOut(e, "GetSeqsCommand", "readList");
		exit(1);
	}
}
//**********************************************************************************************************************
void GetSeqsCommand::readName(){
	try {
	
		string outputFileName = getRootName(namefile) + "pick" +  getExtension(namefile);
		ofstream out;
		openOutputFile(outputFileName, out);

		ifstream in;
		openInputFile(namefile, in);
		string name, firstCol, secondCol;
		
		bool wroteSomething = false;
		
		
		while(!in.eof()){

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
			
			vector<string> validSecond;
			for (int i = 0; i < parsedNames.size(); i++) {
				if (names.count(parsedNames[i]) == 1) {
					validSecond.push_back(parsedNames[i]);
				}
			}

			
			//if the name in the first column is in the set then print it and any other names in second column also in set
			if (names.count(firstCol) == 1) {
			
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
			
			gobble(in);
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) {
			mothurOut("Your file does not contain any sequence from the .accnos file."); mothurOutEndLine();
			remove(outputFileName.c_str()); 
		}
		
	}
	catch(exception& e) {
		errorOut(e, "GetSeqsCommand", "readName");
		exit(1);
	}
}

//**********************************************************************************************************************
void GetSeqsCommand::readGroup(){
	try {
	
		string outputFileName = getRootName(groupfile) + "pick" + getExtension(groupfile);
		ofstream out;
		openOutputFile(outputFileName, out);

		ifstream in;
		openInputFile(groupfile, in);
		string name, group;
		
		bool wroteSomething = false;
		
		while(!in.eof()){

			in >> name;				//read from first column
			in >> group;			//read from second column
			
			//if this name is in the accnos file
			if (names.count(name) == 1) {
				wroteSomething = true;
				
				out << name << '\t' << group << endl;
				
				names.erase(name);
			}
					
			gobble(in);
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) {
			mothurOut("Your file does not contain any sequence from the .accnos file."); mothurOutEndLine();
			remove(outputFileName.c_str()); 
		}

	}
	catch(exception& e) {
		errorOut(e, "GetSeqsCommand", "readGroup");
		exit(1);
	}
}

//**********************************************************************************************************************
//alignreport file has a column header line then all other lines contain 16 columns.  we just want the first column since that contains the name
void GetSeqsCommand::readAlign(){
	try {
		string outputFileName = getRootName(getRootName(alignfile)) + "pick.align.report";
		ofstream out;
		openOutputFile(outputFileName, out);

		ifstream in;
		openInputFile(alignfile, in);
		string name, junk;
		
		bool wroteSomething = false;
		
		//read column headers
		for (int i = 0; i < 16; i++) {  
			if (!in.eof())	{	in >> junk;	 out << junk << '\t';	}
			else			{	break;			}
		}
		out << endl;
		
		while(!in.eof()){

			in >> name;				//read from first column
			
			//if this name is in the accnos file
			if (names.count(name) == 1) {
				wroteSomething = true;
				
				out << name << '\t';
				
				//read rest
				for (int i = 0; i < 15; i++) {  
					if (!in.eof())	{	in >> junk;	 out << junk << '\t';	}
					else			{	break;			}
				}
				out << endl;
				
				names.erase(name);
				
			}else {//still read just don't do anything with it
				//read rest
				for (int i = 0; i < 15; i++) {  
					if (!in.eof())	{	in >> junk;		}
					else			{	break;			}
				}
			}
			
			gobble(in);
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) {
			mothurOut("Your file does not contain any sequence from the .accnos file."); mothurOutEndLine();
			remove(outputFileName.c_str()); 
		}
		
	}
	catch(exception& e) {
		errorOut(e, "GetSeqsCommand", "readAlign");
		exit(1);
	}
}
//**********************************************************************************************************************

void GetSeqsCommand::readAccnos(){
	try {
		
		ifstream in;
		openInputFile(accnosfile, in);
		string name;
		
		while(!in.eof()){
			in >> name;
						
			names.insert(name);
			
			gobble(in);
		}
		in.close();		

	}
	catch(exception& e) {
		errorOut(e, "GetSeqsCommand", "readAccnos");
		exit(1);
	}
}

//**********************************************************************************************************************

