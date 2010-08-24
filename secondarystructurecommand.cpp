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

AlignCheckCommand::AlignCheckCommand(string option)  {
	try {
		abort = false;
		haderror = 0;
			
		//allow user to run help
		if(option == "help") { help(); abort = true; }
		
		else {
			//valid paramters for this command
			string Array[] =  {"fasta","map", "outputdir","inputdir"};
			vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("map");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["map"] = inputDir + it->second;		}
				}
			}

			//check for required parameters
			mapfile = validParameter.validFile(parameters, "map", true);
			if (mapfile == "not open") { abort = true; }
			else if (mapfile == "not found") {  mapfile = "";  m->mothurOut("You must provide an map file."); m->mothurOutEndLine(); abort = true; }	
			
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") {  fastafile = "";  m->mothurOut("You must provide an fasta file."); m->mothurOutEndLine(); abort = true;  }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(fastafile); //if user entered a file with a path then preserve it	
			}

		}

	}
	catch(exception& e) {
		m->errorOut(e, "AlignCheckCommand", "RemoveSeqsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

void AlignCheckCommand::help(){
	try {
		m->mothurOut("The align.check command reads a fasta file and map file.\n");
		m->mothurOut("It outputs a file containing the secondary structure matches in the .align.check file.\n");
		m->mothurOut("The align.check command parameters are fasta and map, both are required.\n");
		m->mothurOut("The align.check command should be in the following format: align.check(fasta=yourFasta, map=yourMap).\n");
		m->mothurOut("Example align.check(map=silva.ss.map, fasta=amazon.fasta).\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCheckCommand", "help");
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
		m->openInputFile(fastafile, in);
		
		ofstream out;
		string outfile = outputDir + m->getRootName(m->getSimpleName(fastafile)) + "align.check";
		m->openOutputFile(outfile, out);
		
		out << "name" << '\t' << "pound" << '\t' << "dash" << '\t' << "plus" << '\t' << "equal" << '\t';
		out << "loop" << '\t' << "tilde" << '\t' << "total" << endl;

		
		while(!in.eof()){
			if (m->control_pressed) { in.close(); out.close(); remove(outfile.c_str()); return 0; }
			
			Sequence seq(in);  m->gobble(in);
			if (seq.getName() != "") {
				statData data = getStats(seq.getAligned());
				
				if (haderror == 1) { break; }
				
				out << seq.getName() << '\t' << data.pound << '\t' << data.dash << '\t' << data.plus << '\t' << data.equal << '\t';
				out << data.loop << '\t' << data.tilde << '\t' << data.total << endl;
			}
		}

		in.close();
		out.close();
		
		if (m->control_pressed) {  remove(outfile.c_str()); return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: "); m->mothurOutEndLine();
		m->mothurOut(outfile); m->mothurOutEndLine();	
		m->mothurOutEndLine();
		
		return 0;		
	}

	catch(exception& e) {
		m->errorOut(e, "AlignCheckCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
void AlignCheckCommand::readMap(){
	try {
			
		structMap.resize(1, 0);
		ifstream in;
		
		m->openInputFile(mapfile, in);
		
		while(!in.eof()){
			int position;
			in >> position;
			structMap.push_back(position);	
			m->gobble(in);
		}
		in.close();

		seqLength = structMap.size();
		
		
		//check you make sure is structMap[10] = 380 then structMap[380] = 10.
		for(int i=0;i<seqLength;i++){
			if(structMap[i] != 0){
				if(structMap[structMap[i]] != i){
					m->mothurOut("Your map file contains an error:  line " + toString(i) + " does not match line " + toString(structMap[i]) + "."); m->mothurOutEndLine();
				}
			}
		}
		
		
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCheckCommand", "readMap");
		exit(1);
	}
}
/**************************************************************************************************/

statData AlignCheckCommand::getStats(string sequence){
	try {
	
		statData data;
		sequence = "*" + sequence; // need to pad the sequence so we can index it by 1
		
		int length = sequence.length();
		
		if (length != seqLength) { m->mothurOut("your sequences are " + toString(length) + " long, but your map file only contains " + toString(seqLength) + " entries. please correct."); m->mothurOutEndLine(); haderror = 1; return data;  }
		
		for(int i=1;i<length;i++){
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
		m->errorOut(e, "AlignCheckCommand", "getStats");
		exit(1);
	}
}


//**********************************************************************************************************************
