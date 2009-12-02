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
		mothurOut("The align.check command reads a fasta file and map file.\n");
		mothurOut("It outputs a file containing the secondary structure matches in the .align.check file.\n");
		mothurOut("The align.check command parameters are fasta and map, both are required.\n");
		mothurOut("The align.check command should be in the following format: align.check(fasta=yourFasta, map=yourMap).\n");
		mothurOut("Example align.check(map=silva.ss.map, fasta=amazon.fasta).\n");
		mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n\n");
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
		
		ifstream in;
		openInputFile(fastafile, in);
		
		ofstream out;
		string outfile = getRootName(fastafile) + "align.check";
		openOutputFile(outfile, out);
		
		out << "name" << '\t' << "pound" << '\t' << "dash" << '\t' << "plus" << '\t' << "equal" << '\t';
		out << "loop" << '\t' << "tilde" << '\t' << "total" << endl;

		
		while(!in.eof()){
			
			Sequence seq(in);  gobble(in);
			if (seq.getName() != "") {
				statData data = getStats(seq.getAligned());
				
				out << seq.getName() << '\t' << data.pound << '\t' << data.dash << '\t' << data.plus << '\t' << data.equal << '\t';
				out << data.loop << '\t' << data.tilde << '\t' << data.total << endl;
			}
		}

		in.close();
		out.close();
		
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
		errorOut(e, "AlignCheckCommand", "readMap");
		exit(1);
	}
}
/**************************************************************************************************/

statData AlignCheckCommand::getStats(string sequence){
	try {
	
		statData data;
		sequence = "*" + sequence; // need to pad the sequence so we can index it by 1
		
		int seqLength = sequence.length();
		for(int i=1;i<seqLength;i++){
			if(structMap[i] != 0){
				if(sequence[i] == 'A'){
					if(sequence[structMap[i]] == 'T')		{	data.tilde++;	}
					else if(sequence[structMap[i]] == 'A')	{	data.pound++;	}
					else if(sequence[structMap[i]] == 'G')	{	data.equal++;	}
					else if(sequence[structMap[i]] == 'C')	{	data.pound++;	}
					else if(sequence[structMap[i]] == '-')	{	data.pound++;	}
					data.total++;
				}
				else if(sequence[i] == 'T'){
					if(sequence[structMap[i]] == 'T')		{	data.plus++;	}
					else if(sequence[structMap[i]] == 'A')	{	data.tilde++;	}
					else if(sequence[structMap[i]] == 'G')	{	data.dash++;	}
					else if(sequence[structMap[i]] == 'C')	{	data.pound++;	}
					else if(sequence[structMap[i]] == '-')	{	data.pound++;	}
					data.total++;
				}
				else if(sequence[i] == 'G'){
					if(sequence[structMap[i]] == 'T')		{	data.dash++;	}
					else if(sequence[structMap[i]] == 'A')	{	data.equal++;	}
					else if(sequence[structMap[i]] == 'G')	{	data.pound++;	}
					else if(sequence[structMap[i]] == 'C')	{	data.tilde++;	}
					else if(sequence[structMap[i]] == '-')	{	data.pound++;	}
					data.total++;
				}
				else if(sequence[i] == 'C'){
					if(sequence[structMap[i]] == 'T')		{	data.pound++;	}
					else if(sequence[structMap[i]] == 'A')	{	data.pound++;	}
					else if(sequence[structMap[i]] == 'G')	{	data.tilde++;	}
					else if(sequence[structMap[i]] == 'C')	{	data.pound++;	}
					else if(sequence[structMap[i]] == '-')	{	data.pound++;	}
					data.total++;
				}
				else if(sequence[i] == '-'){
					if(sequence[structMap[i]] == 'T')		{	data.pound++;	data.total++;	}
					else if(sequence[structMap[i]] == 'A')	{	data.pound++;	data.total++;	}
					else if(sequence[structMap[i]] == 'G')	{	data.pound++;	data.total++;	}
					else if(sequence[structMap[i]] == 'C')	{	data.pound++;	data.total++;	}
					else if(sequence[structMap[i]] == '-')	{		/*donothing*/				}
				}			
			}
			else if(isalnum(sequence[i])){
				data.loop++;
				data.total++;
			}
		}
		return data;
		
	}
	catch(exception& e) {
		errorOut(e, "AlignCheckCommand", "getStats");
		exit(1);
	}
}


//**********************************************************************************************************************
