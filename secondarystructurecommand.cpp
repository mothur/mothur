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
#include "counttable.h"

//**********************************************************************************************************************
vector<string> AlignCheckCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","aligncheck",false,true,true); parameters.push_back(pfasta);
		CommandParameter pmap("map", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(pmap);
		CommandParameter pname("name", "InputTypes", "", "", "namecount", "none", "none","",false,false); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "namecount", "none", "none","",false,false); parameters.push_back(pcount);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCheckCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string AlignCheckCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The align.check command reads a fasta file and map file as well as an optional name or count file.\n";
		helpString += "It outputs a file containing the secondary structure matches in the .align.check file.\n";
		helpString += "The align.check command parameters are fasta and map, both are required.\n";
		helpString += "The align.check command should be in the following format: align.check(fasta=yourFasta, map=yourMap).\n";
		helpString += "Example align.check(map=silva.ss.map, fasta=amazon.fasta).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCheckCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string AlignCheckCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "aligncheck") {  pattern = "[filename],align.check"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "AlignCheckCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
AlignCheckCommand::AlignCheckCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["aligncheck"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "AlignCheckCommand", "AlignCheckCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

AlignCheckCommand::AlignCheckCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		haderror = 0;
			
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
			outputTypes["aligncheck"] = tempOutNames;
			
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
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
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
			mapfile = validParameter.validFile(parameters, "map", true);
			if (mapfile == "not open") { abort = true; }
			else if (mapfile == "not found") {  mapfile = "";  m->mothurOut("You must provide an map file."); m->mothurOutEndLine(); abort = true; }	
			
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not open") { fastafile = ""; abort = true; }
			else if (fastafile == "not found") {  				
				fastafile = m->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { m->setFastaFile(fastafile); }	
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") { namefile = "";  }	
			else { m->setNameFile(namefile); }
			
            countfile = validParameter.validFile(parameters, "count", true);
			if (countfile == "not open") { abort = true; countfile = ""; }	
			else if (countfile == "not found") { countfile = ""; }
			else { m->setCountTableFile(countfile); }
			
            if ((countfile != "") && (namefile != "")) { m->mothurOut("You must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }
            
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(fastafile); //if user entered a file with a path then preserve it	
			}
			
            if (countfile == "") {
                if ((namefile == "") && (fastafile != "")){
                    vector<string> files; files.push_back(fastafile); 
                    parser.getNameFile(files);
                }
            }
		}

	}
	catch(exception& e) {
		m->errorOut(e, "AlignCheckCommand", "AlignCheckCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int AlignCheckCommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		//get secondary structure info.
		readMap();
		
		if (namefile != "") { nameMap = m->readNames(namefile); }
        else if (countfile != "") {
            CountTable ct;
            ct.readTable(countfile, false);
            nameMap = ct.getNameMap();
        }
		
		if (m->control_pressed) { return 0; }
		
		ifstream in;
		m->openInputFile(fastafile, in);
		
		ofstream out;
        map<string, string> variables; 
		variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastafile));
		string outfile = getOutputFileName("aligncheck",variables);
		m->openOutputFile(outfile, out);
		
				
		out << "name" << '\t' << "pound" << '\t' << "dash" << '\t' << "plus" << '\t' << "equal" << '\t';
		out << "loop" << '\t' << "tilde" << '\t' << "total"  << '\t' << "numseqs" << endl;

		vector<int> pound;
		vector<int> dash;
		vector<int> plus;
		vector<int> equal;
		vector<int> loop;
		vector<int> tilde;
		vector<int> total;
		
		int count = 0;
		while(!in.eof()){
			if (m->control_pressed) { in.close(); out.close(); m->mothurRemove(outfile); return 0; }
			
			Sequence seq(in);  m->gobble(in);
			if (seq.getName() != "") {
				statData data = getStats(seq.getAligned());
				
				if (haderror == 1) { m->control_pressed = true; break; }
				
				int num = 1;
				if ((namefile != "") || (countfile != "")) {
					//make sure this sequence is in the namefile, else error 
					map<string, int>::iterator it = nameMap.find(seq.getName());
					
					if (it == nameMap.end()) { m->mothurOut("[ERROR]: " + seq.getName() + " is not in your namefile, please correct."); m->mothurOutEndLine(); m->control_pressed = true; }
					else { num = it->second; }
				}
				
				//for each sequence this sequence represents
				for (int i = 0; i < num; i++) {
					pound.push_back(data.pound);
					dash.push_back(data.dash);
					plus.push_back(data.plus);
					equal.push_back(data.equal);
					loop.push_back(data.loop);
					tilde.push_back(data.tilde);
					total.push_back(data.total);
				}	
				count++;
				
				out << seq.getName() << '\t' << data.pound << '\t' << data.dash << '\t' << data.plus << '\t' << data.equal << '\t';
				out << data.loop << '\t' << data.tilde << '\t' << data.total << '\t' << num << endl;
			}
		}

		in.close();
		out.close();
		
		if (m->control_pressed) {  m->mothurRemove(outfile); return 0; }
		
		sort(pound.begin(), pound.end());
		sort(dash.begin(), dash.end());
		sort(plus.begin(), plus.end());
		sort(equal.begin(), equal.end());
		sort(loop.begin(), loop.end());
		sort(tilde.begin(), tilde.end());
		sort(total.begin(), total.end());
		int size = pound.size();
		
		int ptile0_25	= int(size * 0.025);
		int ptile25		= int(size * 0.250);
		int ptile50		= int(size * 0.500);
		int ptile75		= int(size * 0.750);
		int ptile97_5	= int(size * 0.975);
		int ptile100	= size - 1;
		
		if (m->control_pressed) {  m->mothurRemove(outfile); return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("\t\tPound\tDash\tPlus\tEqual\tLoop\tTilde\tTotal"); m->mothurOutEndLine();
		m->mothurOut("Minimum:\t" + toString(pound[0]) + "\t" + toString(dash[0]) + "\t" + toString(plus[0]) + "\t" + toString(equal[0]) + "\t" + toString(loop[0]) + "\t" + toString(tilde[0]) + "\t" + toString(total[0])); m->mothurOutEndLine();
		m->mothurOut("2.5%-tile:\t" + toString(pound[ptile0_25]) + "\t" + toString(dash[ptile0_25]) + "\t" + toString(plus[ptile0_25]) + "\t" + toString(equal[ptile0_25]) + "\t"+ toString(loop[ptile0_25]) + "\t"+ toString(tilde[ptile0_25]) + "\t"+ toString(total[ptile0_25])); m->mothurOutEndLine();
		m->mothurOut("25%-tile:\t" + toString(pound[ptile25]) + "\t" + toString(dash[ptile25]) + "\t" + toString(plus[ptile25]) + "\t" + toString(equal[ptile25]) + "\t" + toString(loop[ptile25]) + "\t" + toString(tilde[ptile25]) + "\t" + toString(total[ptile25])); m->mothurOutEndLine();
		m->mothurOut("Median: \t" + toString(pound[ptile50]) + "\t" + toString(dash[ptile50]) + "\t" + toString(plus[ptile50]) + "\t" + toString(equal[ptile50]) + "\t" + toString(loop[ptile50]) + "\t" + toString(tilde[ptile50]) + "\t" + toString(total[ptile50])); m->mothurOutEndLine();
		m->mothurOut("75%-tile:\t" + toString(pound[ptile75]) + "\t" + toString(dash[ptile75]) + "\t" + toString(plus[ptile75]) + "\t" + toString(equal[ptile75]) + "\t" + toString(loop[ptile75]) + "\t" + toString(tilde[ptile75]) + "\t" + toString(total[ptile75])); m->mothurOutEndLine();
		m->mothurOut("97.5%-tile:\t" + toString(pound[ptile97_5]) + "\t" + toString(dash[ptile97_5]) + "\t" + toString(plus[ptile97_5]) + "\t" + toString(equal[ptile97_5]) + "\t" + toString(loop[ptile97_5]) + "\t" + toString(tilde[ptile97_5]) + "\t" + toString(total[ptile97_5])); m->mothurOutEndLine();
		m->mothurOut("Maximum:\t" + toString(pound[ptile100]) + "\t" + toString(dash[ptile100]) + "\t" + toString(plus[ptile100]) + "\t" + toString(equal[ptile100]) + "\t" + toString(loop[ptile100]) + "\t" + toString(tilde[ptile100]) + "\t" + toString(total[ptile100])); m->mothurOutEndLine();
		if ((namefile == "") && (countfile == "")) {  m->mothurOut("# of Seqs:\t" + toString(count)); m->mothurOutEndLine(); }
		else { m->mothurOut("# of unique seqs:\t" + toString(count)); m->mothurOutEndLine(); m->mothurOut("total # of seqs:\t" + toString(size)); m->mothurOutEndLine(); }
		
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		m->mothurOut(outfile); m->mothurOutEndLine();	outputNames.push_back(outfile); outputTypes["aligncheck"].push_back(outfile);
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
