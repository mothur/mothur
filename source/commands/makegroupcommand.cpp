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
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
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
				if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["group"] = tempOutNames;
		
            //if the user changes the output directory command factory will send this info to us in the output parameter
            outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";		}
            
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}

			fastaFileName = validParameter.valid(parameters, "fasta");
			if (fastaFileName == "not found") { 				//if there is a current fasta file, use it
				string filename = current->getFastaFile(); 
				if (filename != "") { fastaFileNames.push_back(filename); m->mothurOut("Using " + filename + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { 
				util.splitAtDash(fastaFileName, fastaFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < fastaFileNames.size(); i++) {
					
					bool ignore = false;
					if (fastaFileNames[i] == "current") { 
						fastaFileNames[i] = current->getFastaFile(); 
						if (fastaFileNames[i] != "") {  
							m->mothurOut("Using " + fastaFileNames[i] + " as input file for the fasta parameter where you had given current."); m->mothurOutEndLine(); 
							filename += util.getRootName(util.getSimpleName(fastaFileNames[i]));
						}
						else { 	
							m->mothurOut("You have no current fastafile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							fastaFileNames.erase(fastaFileNames.begin()+i);
							i--;
						}
					}
					
                    if (!ignore) {
                        if (util.checkLocations(fastaFileNames[i], current->getLocations())) { current->setFastaFile(fastaFileNames[i]); }
                        else { fastaFileNames.erase(fastaFileNames.begin()+i); i--; } //erase from file list
                    }
				}
				
				//prevent giantic file name
                map<string, string> variables; 
                variables["[filename]"] = filename;
				if (fastaFileNames.size() > 3) { variables["[filename]"] = "merge."; }
				filename = getOutputFileName("group",variables);  
				
				//make sure there is at least one valid file left
				if (fastaFileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
			
			output = validParameter.valid(parameters, "output");			
			if (output == "not found") { output = "";  }
			else{ filename = output; }
			
			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { m->mothurOut("groups is a required parameter for the make.group command."); m->mothurOutEndLine(); abort = true;  }
			else { util.splitAtDash(groups, groupsNames);	}

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
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		if (outputDir == "") { outputDir = util.hasPath(fastaFileNames[0]); }
			
		filename = outputDir + filename;
		
		ofstream out;
		util.openOutputFile(filename, out);
		
		for (int i = 0; i < fastaFileNames.size(); i++) {
		
			if (m->getControl_pressed()) { outputTypes.clear(); out.close(); util.mothurRemove(filename); return 0; }
			
			ifstream in;
			util.openInputFile(fastaFileNames[i], in);
			
			while (!in.eof()) {
				
				Sequence seq(in); util.gobble(in);
				
				if (m->getControl_pressed()) { outputTypes.clear();  in.close(); out.close(); util.mothurRemove(filename); return 0; }
				
				if (seq.getName() != "") {	out << seq.getName() << '\t' << groupsNames[i] << endl;		}
			}
			in.close();
		}
		
		out.close();
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: " + filename); m->mothurOutEndLine(); outputNames.push_back(filename); outputTypes["group"].push_back(filename); 
		m->mothurOutEndLine();
		
		//set group file as new current groupfile
		string currentName = "";
		itTypes = outputTypes.find("group");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setGroupFile(currentName); }
		}
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeGroupCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************


