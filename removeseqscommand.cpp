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

//**********************************************************************************************************************

RemoveSeqsCommand::RemoveSeqsCommand(string option){
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta","name", "group", "alignreport", "accnos" };
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
			
			if ((fastafile == "") && (namefile == "") && (groupfile == "") && (alignfile == ""))  { mothurOut("You must provide one of the following: fasta, name, group, alignreport."); mothurOutEndLine(); abort = true; }
			
			if (parameters.size() > 2) { mothurOut("You may only enter one of the following: fasta, name, group, alignreport."); mothurOutEndLine(); abort = true;  }
		}

	}
	catch(exception& e) {
		errorOut(e, "RemoveSeqsCommand", "RemoveSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void RemoveSeqsCommand::help(){
	try {
		mothurOut("The remove.seqs command reads an .accnos file and one of the following file types: fasta, name, group or alignreport file.\n");
		mothurOut("It outputs a file containing the sequences NOT in the .accnos file.\n");
		mothurOut("The remove.seqs command parameters are accnos, fasta, name, group and alignreport.  You must provide accnos and one of the other parameters.\n");
		mothurOut("The remove.seqs command should be in the following format: remove.seqs(accnos=yourAccnos, fasta=yourFasta).\n");
		mothurOut("Example remove.seqs(accnos=amazon.accnos, fasta=amazon.fasta).\n");
		mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n\n");
	}
	catch(exception& e) {
		errorOut(e, "RemoveSeqsCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

int RemoveSeqsCommand::execute(){
	try {
		
		if (abort == true) { return 0; }
		
		//get names you want to keep
		readAccnos();
		
		//read through the correct file and output lines you want to keep
		if (fastafile != "")		{		readFasta();	}
		else if (namefile != "")	{		readName();		}
		else if (groupfile != "")	{		readGroup();	}
		else if (alignfile != "")	{		readAlign();	}
		
		return 0;		
	}

	catch(exception& e) {
		errorOut(e, "RemoveSeqsCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
void RemoveSeqsCommand::readFasta(){
	try {
		string outputFileName = getRootName(fastafile) + "pick" + getExtension(fastafile);
		ofstream out;
		openOutputFile(outputFileName, out);
		
		ifstream in;
		openInputFile(fastafile, in);
		string name;
		
		bool wroteSomething = false;
		
		while(!in.eof()){
			Sequence currSeq(in);
			name = currSeq.getName();
			
			//if this name is in the accnos file
			if (names.count(name) == 0) {
				wroteSomething = true;
				
				currSeq.printSequence(out);
			}else {		names.erase(name);		}
			
			gobble(in);
		}
		in.close();	
		out.close();
		
		if (wroteSomething == false) {
			mothurOut("Your file contains only sequences from the .accnos file."); mothurOutEndLine();
			remove(outputFileName.c_str()); 
		}

	}
	catch(exception& e) {
		errorOut(e, "RemoveSeqsCommand", "readFasta");
		exit(1);
	}
}

//**********************************************************************************************************************
void RemoveSeqsCommand::readName(){
	try {
	
		string outputFileName = getRootName(namefile) + "pick" + getExtension(namefile);

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

			vector<string> validSecond;  validSecond.clear();
			for (int i = 0; i < parsedNames.size(); i++) {
				if (names.count(parsedNames[i]) == 0) {
					validSecond.push_back(parsedNames[i]);
				}
			}
			
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
			
			gobble(in);
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) {
			mothurOut("Your file contains only sequences from the .accnos file."); mothurOutEndLine();
			remove(outputFileName.c_str()); 
		}
		
	}
	catch(exception& e) {
		errorOut(e, "RemoveSeqsCommand", "readName");
		exit(1);
	}
}

//**********************************************************************************************************************
void RemoveSeqsCommand::readGroup(){
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
			if (names.count(name) == 0) {
				wroteSomething = true;
				out << name << '\t' << group << endl;
			}else {		names.erase(name);		}
					
			gobble(in);
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) {
			mothurOut("Your file contains only sequences from the .accnos file."); mothurOutEndLine();
			remove(outputFileName.c_str()); 
		}

	}
	catch(exception& e) {
		errorOut(e, "RemoveSeqsCommand", "readGroup");
		exit(1);
	}
}

//**********************************************************************************************************************
//alignreport file has a column header line then all other lines contain 16 columns.  we just want the first column since that contains the name
void RemoveSeqsCommand::readAlign(){
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
				names.erase(name);	
				
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
			mothurOut("Your file contains only sequences from the .accnos file."); mothurOutEndLine();
			remove(outputFileName.c_str()); 
		}
		
	}
	catch(exception& e) {
		errorOut(e, "RemoveSeqsCommand", "readAlign");
		exit(1);
	}
}
//**********************************************************************************************************************
void RemoveSeqsCommand::readAccnos(){
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
		errorOut(e, "RemoveSeqsCommand", "readAccnos");
		exit(1);
	}
}

//**********************************************************************************************************************


