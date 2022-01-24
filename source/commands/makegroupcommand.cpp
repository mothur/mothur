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
		
        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["group"] = tempOutNames;
        
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
MakeGroupCommand::MakeGroupCommand(string option) : Command()  {
	try {
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string, string> parameters = parser.getParameters(); 
			
			ValidParameters validParameter;
            
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.validPath(parameters, "inputdir");
			if (inputDir == "not found"){	inputDir = "";		}

            fastaFileNames = validParameter.validFiles(parameters, "fasta");
            if (fastaFileNames.size() != 0) {
                if (fastaFileNames[0] == "not open") { abort = true; }
                else {
                    current->setFastaFile(fastaFileNames[0]);
                    filename = util.getRootName(util.getSimpleName(fastaFileNames[0]));
                    //prevent giantic file name
                    map<string, string> variables;
                    variables["[filename]"] = filename;
                    if (fastaFileNames.size() > 3) { variables["[filename]"] = "merge."; }
                    filename = getOutputFileName("group",variables);
                }
            }
            
            //make sure there is at least one valid file left
            if (fastaFileNames.size() == 0) { m->mothurOut("[ERROR]: no valid files.\n");  abort = true; }
            
			output = validParameter.validPath(parameters, "output");
			if (output == "not found") { output = "";  }
			else{ filename = output; }
			
			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { m->mothurOut("groups is a required parameter for the make.group command.\n");  abort = true;  }
			else { util.splitAtDash(groups, groupsNames);	}

			if (groupsNames.size() != fastaFileNames.size()) { m->mothurOut("You do not have the same number of valid fastfile files as groups.  This could be because we could not open a fastafile.\n");  abort = true;  }
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
		
		if (outputdir == "") { outputdir = util.hasPath(fastaFileNames[0]); }
			
		filename = outputdir + filename;
		
		ofstream out;
		util.openOutputFile(filename, out);
		
		for (int i = 0; i < fastaFileNames.size(); i++) {
		
			if (m->getControl_pressed()) { outputTypes.clear(); out.close(); util.mothurRemove(filename); return 0; }
			
			ifstream in; util.openInputFile(fastaFileNames[i], in);
			
			while (!in.eof()) {
				
				Sequence seq(in); util.gobble(in);
				
				if (m->getControl_pressed()) { outputTypes.clear();  in.close(); out.close(); util.mothurRemove(filename); return 0; }
				
				if (seq.getName() != "") {	out << seq.getName() << '\t' << groupsNames[i] << endl;		}
			}
			in.close();
		}
		
		out.close();
		
		
		m->mothurOut("\nOutput File Names: " + filename + "\n\n");
        outputNames.push_back(filename); outputTypes["group"].push_back(filename);
		
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


