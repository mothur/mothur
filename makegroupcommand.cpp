/*
 *  makegroupcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 5/7/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "makegroupcommand.h"
#include "sequence.hpp"

//**********************************************************************************************************************
vector<string> MakeGroupCommand::getValidParameters(){	
	try {
		string Array[] =  {"fasta", "output","groups","outputdir","inputdir"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeGroupCommand", "getValidParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
MakeGroupCommand::MakeGroupCommand(){	
	try {
		abort = true; calledHelp = true; 
		vector<string> tempOutNames;
		outputTypes["group"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeGroupCommand", "MakeGroupCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> MakeGroupCommand::getRequiredParameters(){	
	try {
		string Array[] =  {"fasta","groups"};
		vector<string> myArray (Array, Array+(sizeof(Array)/sizeof(string)));
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeGroupCommand", "getRequiredParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> MakeGroupCommand::getRequiredFiles(){	
	try {
		vector<string> myArray;
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeGroupCommand", "getRequiredFiles");
		exit(1);
	}
}
//**********************************************************************************************************************

MakeGroupCommand::MakeGroupCommand(string option)  {
	try {
		
		abort = false; calledHelp = false;   
	
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		
		else {
			
			//valid paramters for this command
			string AlignArray[] =  {"fasta","groups","output","outputdir","inputdir"};
			vector<string> myArray (AlignArray, AlignArray+(sizeof(AlignArray)/sizeof(string)));
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters(); 
			
			ValidParameters validParameter;
			map<string, string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["group"] = tempOutNames;
		
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}

			fastaFileName = validParameter.validFile(parameters, "fasta", false);
			if (fastaFileName == "not found") { m->mothurOut("fasta is a required parameter for the make.group command."); m->mothurOutEndLine(); abort = true;  }
			else { 
				m->splitAtDash(fastaFileName, fastaFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < fastaFileNames.size(); i++) {
					if (inputDir != "") {
						string path = m->hasPath(fastaFileNames[i]);
						//if the user has not given a path then, add inputdir. else leave path alone.
						if (path == "") {	fastaFileNames[i] = inputDir + fastaFileNames[i];		}
					}
	
					ifstream in;
					int ableToOpen = m->openInputFile(fastaFileNames[i], in, "noerror");
				
					//if you can't open it, try default location
					if (ableToOpen == 1) {
						if (m->getDefaultPath() != "") { //default path is set
							string tryPath = m->getDefaultPath() + m->getSimpleName(fastaFileNames[i]);
							m->mothurOut("Unable to open " + fastaFileNames[i] + ". Trying default " + tryPath); m->mothurOutEndLine();
							ifstream in2;
							ableToOpen = m->openInputFile(tryPath, in2, "noerror");
							in2.close();
							fastaFileNames[i] = tryPath;
						}
					}
					
					//if you can't open it, try default location
					if (ableToOpen == 1) {
						if (m->getOutputDir() != "") { //default path is set
							string tryPath = m->getOutputDir() + m->getSimpleName(fastaFileNames[i]);
							m->mothurOut("Unable to open " + fastaFileNames[i] + ". Trying output directory " + tryPath); m->mothurOutEndLine();
							ifstream in2;
							ableToOpen = m->openInputFile(tryPath, in2, "noerror");
							in2.close();
							fastaFileNames[i] = tryPath;
						}
					}
					in.close();
					
					if (ableToOpen == 1) { 
						m->mothurOut("Unable to open " + fastaFileNames[i] + ". It will be disregarded."); m->mothurOutEndLine();
						//erase from file list
						fastaFileNames.erase(fastaFileNames.begin()+i);
						i--;
					}else{  filename += m->getRootName(m->getSimpleName(fastaFileNames[i]));  }
				}
				
				//prevent giantic file name
				if (fastaFileNames.size() > 3) { filename = outputDir + "merge.groups"; }
				else {  filename += "groups";  }
				
				//make sure there is at least one valid file left
				if (fastaFileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
			output = validParameter.validFile(parameters, "output", false);			
			if (output == "not found") { output = "";  }
			else{ filename = output; }
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { m->mothurOut("groups is a required parameter for the make.group command."); m->mothurOutEndLine(); abort = true;  }
			else { m->splitAtDash(groups, groupsNames);	}

			if (groupsNames.size() != fastaFileNames.size()) { m->mothurOut("You do not have the same number of valid fastfile files as groups.  This could be because we could not open a fastafile."); m->mothurOutEndLine(); abort = true;  }
		}
		cout << "here" << endl;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeGroupCommand", "MakeGroupCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

MakeGroupCommand::~MakeGroupCommand(){	}

//**********************************************************************************************************************

void MakeGroupCommand::help(){
	try {
		m->mothurOut("The make.group command reads a fasta file or series of fasta files and creates a groupfile.\n");
		m->mothurOut("The make.group command parameters are fasta, groups and output. Fasta and group are required.\n");
		m->mothurOut("The output parameter allows you to specify the name of groupfile created. \n");
		m->mothurOut("The make.group command should be in the following format: \n");
		m->mothurOut("make.group(fasta=yourFastaFiles, groups=yourGroups). \n");
		m->mothurOut("Example make.group(fasta=seqs1.fasta-seq2.fasta-seqs3.fasta, groups=A-B-C)\n");
		m->mothurOut("Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFiles).\n\n");
	}
	catch(exception& e) {
		m->errorOut(e, "MakeGroupCommand", "help");
		exit(1);
	}
}


//**********************************************************************************************************************

int MakeGroupCommand::execute(){
	try {
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		if (outputDir == "") { outputDir = m->hasPath(fastaFileNames[0]); }
			
		filename = outputDir + filename;
		
		ofstream out;
		m->openOutputFile(filename, out);
		
		for (int i = 0; i < fastaFileNames.size(); i++) {
		
			if (m->control_pressed) { outputTypes.clear(); out.close(); remove(filename.c_str()); return 0; }
			
			ifstream in;
			m->openInputFile(fastaFileNames[i], in);
			
			while (!in.eof()) {
				
				Sequence seq(in, "no align"); m->gobble(in);
				
				if (m->control_pressed) { outputTypes.clear();  in.close(); out.close(); remove(filename.c_str()); return 0; }
				
				if (seq.getName() != "") {	out << seq.getName() << '\t' << groupsNames[i] << endl;		}
			}
			in.close();
		}
		
		out.close();
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Name: " + filename); m->mothurOutEndLine(); outputNames.push_back(filename); outputTypes["group"].push_back(filename); 
		m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeGroupCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************


