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
vector<string> MakeGroupCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","group",false,true,true); parameters.push_back(pfasta);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false,true); parameters.push_back(pgroups);
		CommandParameter poutput("output", "String", "", "", "", "", "","",false,false); parameters.push_back(poutput);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeGroupCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string MakeGroupCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The make.group command reads a fasta file or series of fasta files and creates a groupfile.\n";
		helpString += "The make.group command parameters are fasta, groups and output. Fasta and group are required.\n";
		helpString += "The output parameter allows you to specify the name of groupfile created. \n";
		helpString += "The make.group command should be in the following format: \n";
		helpString += "make.group(fasta=yourFastaFiles, groups=yourGroups). \n";
		helpString += "Example make.group(fasta=seqs1.fasta-seq2.fasta-seqs3.fasta, groups=A-B-C)\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFiles).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeGroupCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string MakeGroupCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "group") {  pattern = "[filename],groups"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeGroupCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
MakeGroupCommand::MakeGroupCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["group"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeGroupCommand", "MakeGroupCommand");
		exit(1);
	}
}

//**********************************************************************************************************************

MakeGroupCommand::MakeGroupCommand(string option)  {
	try {
		
		abort = false; calledHelp = false;   
	
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
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
		
            //if the user changes the output directory command factory will send this info to us in the output parameter
            outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
            
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}

			fastaFileName = validParameter.validFile(parameters, "fasta", false);
			if (fastaFileName == "not found") { 				//if there is a current fasta file, use it
				string filename = m->getFastaFile(); 
				if (filename != "") { fastaFileNames.push_back(filename); m->mothurOut("Using " + filename + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { 
				m->splitAtDash(fastaFileName, fastaFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < fastaFileNames.size(); i++) {
					
					bool ignore = false;
					if (fastaFileNames[i] == "current") { 
						fastaFileNames[i] = m->getFastaFile(); 
						if (fastaFileNames[i] != "") {  
							m->mothurOut("Using " + fastaFileNames[i] + " as input file for the fasta parameter where you had given current."); m->mothurOutEndLine(); 
							filename += m->getRootName(m->getSimpleName(fastaFileNames[i]));
						}
						else { 	
							m->mothurOut("You have no current fastafile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							fastaFileNames.erase(fastaFileNames.begin()+i);
							i--;
						}
					}
					
					if (!ignore) {
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
						}else{  filename += m->getRootName(m->getSimpleName(fastaFileNames[i]));  m->setFastaFile(fastaFileNames[i]); }
					}
				}
				
				//prevent giantic file name
                map<string, string> variables; 
                variables["[filename]"] = filename;
				if (fastaFileNames.size() > 3) { variables["[filename]"] = outputDir + "merge"; }
				filename = getOutputFileName("group",variables);  
				
				//make sure there is at least one valid file left
				if (fastaFileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
			
			output = validParameter.validFile(parameters, "output", false);			
			if (output == "not found") { output = "";  }
			else{ filename = output; }
			
			groups = validParameter.validFile(parameters, "groups", false);			
			if (groups == "not found") { m->mothurOut("groups is a required parameter for the make.group command."); m->mothurOutEndLine(); abort = true;  }
			else { m->splitAtDash(groups, groupsNames);	}

			if (groupsNames.size() != fastaFileNames.size()) { m->mothurOut("You do not have the same number of valid fastfile files as groups.  This could be because we could not open a fastafile."); m->mothurOutEndLine(); abort = true;  }
		}
	}
	catch(exception& e) {
		m->errorOut(e, "MakeGroupCommand", "MakeGroupCommand");
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
		
			if (m->control_pressed) { outputTypes.clear(); out.close(); m->mothurRemove(filename); return 0; }
			
			ifstream in;
			m->openInputFile(fastaFileNames[i], in);
			
			while (!in.eof()) {
				
				Sequence seq(in, "no align"); m->gobble(in);
				
				if (m->control_pressed) { outputTypes.clear();  in.close(); out.close(); m->mothurRemove(filename); return 0; }
				
				if (seq.getName() != "") {	out << seq.getName() << '\t' << groupsNames[i] << endl;		}
			}
			in.close();
		}
		
		out.close();
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: " + filename); m->mothurOutEndLine(); outputNames.push_back(filename); outputTypes["group"].push_back(filename); 
		m->mothurOutEndLine();
		
		//set group file as new current groupfile
		string current = "";
		itTypes = outputTypes.find("group");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setGroupFile(current); }
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeGroupCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************


