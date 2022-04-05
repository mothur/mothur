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
#include "counttable.h"


//**********************************************************************************************************************
vector<string> MakeGroupCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","group",false,true,true); parameters.push_back(pfasta);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false,true); parameters.push_back(pgroups);
		CommandParameter poutput("output", "String", "", "", "", "", "","",false,false); parameters.push_back(poutput);
        CommandParameter pformat("format", "Multiple", "count-group", "count", "", "", "","",false,false,true); parameters.push_back(pformat);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["group"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
        
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
		helpString += "The make.group (also called make.count) command reads a fasta file or series of fasta files and creates a group file or count file.\n";
		helpString += "The make.group command parameters are fasta, groups, format and output. Fasta and groups are required.\n";
		helpString += "The output parameter allows you to specify the name of group file or count file created. \n";
        helpString += "The format parameter allows you to specify whether the outputtted file is a group file or count file. Default=count. \n";
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
        else if (type == "count") {  pattern = "[filename],count_table"; }
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
                else {  current->setFastaFile(fastaFileNames[0]); }
            }
            
            //make sure there is at least one valid file left
            if (fastaFileNames.size() == 0) { m->mothurOut("[ERROR]: no valid files.\n");  abort = true; }
            
			output = validParameter.validPath(parameters, "output");
			if (output == "not found") { output = "";  }
			
            format = validParameter.valid(parameters, "format");        if (format == "not found"){    format = "count";    }
            if ((format != "count") && (format != "group")) { m->mothurOut("\n[WARNING]: invalid format option: choices are count or group, using count.\n");  format="count"; }
            
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
        
        util.checkGroupNames(groupsNames);
		
        map<string, string> seqGroup; map<string, long long> groupCounts;
		for (int i = 0; i < fastaFileNames.size(); i++) {
            
            if (m->getControl_pressed()) { break; }
            
            m->mothurOutJustToScreen("\nAssigning sequences from file " + fastaFileNames[i] + " to group " +  groupsNames[i] + ":\t");
		
			ifstream in; util.openInputFile(fastaFileNames[i], in);
			
            long long count = 0;
			while (!in.eof()) {
				
                if (m->getControl_pressed()) { break; }
                
				Sequence seq(in); gobble(in);
								
                if (seq.getName() != "") {	seqGroup[seq.getName()] = groupsNames[i];	count++;	}
			}
			in.close();
            m->mothurOutJustToScreen(toString(count) + " sequences assigned to group " +  groupsNames[i] + "\n");
            groupCounts[groupsNames[i]] = count;
		}
        
        if (m->getControl_pressed()) { return 0; }
        
        //if user provided output filename, then use it
        string outputFileName = util.getFullPathName(output);
        
        //if no output filename given, create one
        if (output == "") {
            map<string, string> variables;
            if (outputdir == "") { outputdir = util.hasPath(fastaFileNames[0]); }
            variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(fastaFileNames[0]));
            if (fastaFileNames.size() > 1) { variables["[filename]"] = outputdir + "merge."; }
            outputFileName = getOutputFileName(format,variables);
        }
        outputNames.push_back(outputFileName); outputTypes[format].push_back(outputFileName);
        
        if (format == "count") {
            CountTable ct; ct.createTable(seqGroup);
            ct.printCompressedTable(outputFileName);
        }else{
            ofstream out; util.openOutputFile(outputFileName, out);
            for (map<string, string>::iterator it = seqGroup.begin(); it != seqGroup.end(); it++) {
                out << it->first << '\t' << it->second << endl;
            }
            out.close();
        }
        
        long long total = 0;
        if (groupCounts.size() != 0) {  m->mothurOut("\nGroup count: \n");  }
        for (map<string, long long>::iterator it = groupCounts.begin(); it != groupCounts.end(); it++) { total += it->second; m->mothurOut(it->first + "\t" + toString(it->second) + "\n"); }
        if (total != 0) { m->mothurOut("\nTotal of all groups is " + toString(total) + "\n"); }

		m->mothurOut("\nOutput File Names: " + outputFileName + "\n\n");
        
		//set group file as new current groupfile
		string currentName = "";
		itTypes = outputTypes.find("group");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setGroupFile(currentName); }
		}
        
        itTypes = outputTypes.find("count");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setCountFile(currentName); }
        }
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeGroupCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************


