/*
 *  secondarystructurecommand.cpp
 *  Mothur
 *
 *  Created by westcott on 9/18/09.
 *  Copyright 2009 Schloss Lab. All rights reserved.
 *
 */

#include "secondarystructurecommand.h"
#include "sequence.hpp"

//**********************************************************************************************************************

AlignCheckCommand::AlignCheckCommand(string option){
	try {
		abort = false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta","map"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			for (map<string,string>::iterator it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//check for required parameters
			mapfile = validParameter.validFile(parameters, "map", true);
			if (mapfile == "not open") { abort = true; }
			else if (mapfile == "not found") {  mapfile = "";  mothurOut("You must provide an map file."); mothurOutEndLine(); abort = true; }	
			
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") {  fastafile = "";  mothurOut("You must provide an fasta file."); mothurOutEndLine(); abort = true;  }	
			
		}

	}
	catch(exception& e) {
		errorOut(e, "AlignCheckCommand", "RemoveSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void AlignCheckCommand::help(){
	try {
		//mothurOut("The remove.seqs command reads an .accnos file and one of the following file types: fasta, name, group or alignreport file.\n");
		//mothurOut("It outputs a file containing the sequences NOT in the .accnos file.\n");
		//mothurOut("The remove.seqs command parameters are accnos, fasta, name, group and alignreport.  You must provide accnos and one of the other parameters.\n");
		//mothurOut("The remove.seqs command should be in the following format: remove.seqs(accnos=yourAccnos, fasta=yourFasta).\n");
		//mothurOut("Example remove.seqs(accnos=amazon.accnos, fasta=amazon.fasta).\n");
		//mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n\n");
	}
	catch(exception& e) {
		errorOut(e, "AlignCheckCommand", "help");
		exit(1);
	}
}

//**********************************************************************************************************************

int AlignCheckCommand::execute(){
	try {
		
		if (abort == true) { return 0; }
		
		//get secondary structure info.
		readMap();
		
	
		
		return 0;		
	}

	catch(exception& e) {
		errorOut(e, "AlignCheckCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
void AlignCheckCommand::readMap(){
	try {
			
		structMap.resize(1, 0);
		ifstream in;
		
		openInputFile(mapfile, in);
		
		while(!in.eof()){
			int position;
			in >> position;
			structMap.push_back(position);	
			gobble(in);
		}
		in.close();
		
		seqLength = structMap.size();
		
		
		//check you make sure is structMap[10] = 380 then structMap[380] = 10.
		for(int i=0;i<seqLength;i++){
			if(structMap[i] != 0){
				if(structMap[structMap[i]] != i){
					mothurOut("Your map file contains an error:  line " + toString(i) + " does not match line " + toString(structMap[i]) + "."); mothurOutEndLine();
				}
			}
		}
		
		
	}
	catch(exception& e) {
		errorOut(e, "AlignCheckCommand", "readFasta");
		exit(1);
	}
}

//**********************************************************************************************************************
