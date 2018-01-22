/*
 *  countgroupscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 8/9/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "countgroupscommand.h"

#include "inputdata.h"

//**********************************************************************************************************************
vector<string> CountGroupsCommand::setParameters(){	
	try {
		CommandParameter pshared("shared", "InputTypes", "", "", "sharedGroup", "sharedGroup", "none","summary",false,false,true); parameters.push_back(pshared);
		CommandParameter pgroup("group", "InputTypes", "", "", "sharedGroup", "sharedGroup", "none","summary",false,false,true); parameters.push_back(pgroup);
        CommandParameter pcount("count", "InputTypes", "", "", "sharedGroup", "sharedGroup", "none","summary",false,false,true); parameters.push_back(pcount);
		CommandParameter paccnos("accnos", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(paccnos);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "CountGroupsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string CountGroupsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "summary") {  pattern = "[filename],count.summary"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "PrimerDesignCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
string CountGroupsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The count.groups command counts sequences from a specific group or set of groups from the following file types: group, count or shared file.\n";
		helpString += "The count.groups command parameters are accnos, group, shared and groups. You must provide a group or shared file.\n";
		helpString += "The accnos parameter allows you to provide a file containing the list of groups.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like.  You can separate group names with dashes.\n";
		helpString += "The count.groups command should be in the following format: count.groups(accnos=yourAccnos, group=yourGroupFile).\n";
		helpString += "Example count.groups(accnos=amazon.accnos, group=amazon.groups).\n";
		helpString += "or count.groups(groups=pasture, group=amazon.groups).\n";
		helpString += "Note: No spaces between parameter labels (i.e. group), '=' and parameters (i.e.yourGroupFile).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "CountGroupsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
CountGroupsCommand::CountGroupsCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
        vector<string> tempOutNames;
		outputTypes["summary"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "CountGroupsCommand", "CountGroupsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
CountGroupsCommand::CountGroupsCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
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
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";		}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("accnos");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["accnos"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
				
				it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
			}
			
            vector<string> tempOutNames;
            outputTypes["summary"] = tempOutNames;
			
			//check for required parameters
			accnosfile = validParameter.validFile(parameters, "accnos");
			if (accnosfile == "not open") { abort = true; }
			else if (accnosfile == "not found") {  accnosfile = ""; }
			else { current->setAccnosFile(accnosfile); }
			
			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = ""; }
			else {
				util.splitAtDash(groups, Groups);
                if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
			}
			
			sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }
			else if (sharedfile == "not found") {  sharedfile = "";  }
			else { current->setSharedFile(sharedfile); }
			
			groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not open") { groupfile = ""; abort = true; }
			else if (groupfile == "not found") {  	groupfile = "";	}
			else { current->setGroupFile(groupfile); }
            
            countfile = validParameter.validFile(parameters, "count");
            if (countfile == "not open") { countfile = ""; abort = true; }
            else if (countfile == "not found") { countfile = "";  }	
            else { 
                current->setCountFile(countfile); 
                CountTable ct;
                if (!ct.testGroups(countfile)) { m->mothurOut("[ERROR]: Your count file does not have any group information, aborting."); m->mothurOutEndLine(); abort=true; }
            }
            
            if ((groupfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: group or count."); m->mothurOutEndLine(); abort=true;
            }

			
			if ((sharedfile == "") && (groupfile == "") && (countfile == "")) { 
				//give priority to shared, then group
				sharedfile = current->getSharedFile(); 
				if (sharedfile != "") {  m->mothurOut("Using " + sharedfile + " as input file for the shared parameter."); m->mothurOutEndLine(); }
				else { 
					groupfile = current->getGroupFile(); 
					if (groupfile != "") { m->mothurOut("Using " + groupfile + " as input file for the group parameter."); m->mothurOutEndLine(); }
					else { 
						countfile = current->getCountFile(); 
                        if (countfile != "") { m->mothurOut("Using " + countfile + " as input file for the count parameter."); m->mothurOutEndLine(); }
                        else { 
                            m->mothurOut("You have no current groupfile, countfile or sharedfile and one is required."); m->mothurOutEndLine(); abort = true;
                        }
					}
				}
			}
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "CountGroupsCommand", "CountGroupsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int CountGroupsCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		//get groups you want to remove
		if (accnosfile != "") { util.readAccnos(accnosfile, Groups);  }
		
		if (groupfile != "") {
            map<string, string> variables; 
            string thisOutputDir = outputDir;
            if (outputDir == "") {  thisOutputDir += util.hasPath(groupfile);  }
            variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(groupfile));
            string outputFileName = getOutputFileName("summary", variables);
            outputNames.push_back(outputFileName); outputTypes["summary"].push_back(outputFileName);
            ofstream out;
            util.openOutputFile(outputFileName, out);
            
			GroupMap groupMap(groupfile);
			groupMap.readMap();
			
			vector<string> nameGroups = groupMap.getNamesOfGroups();
            if (Groups.size() == 0) { Groups = nameGroups;  }
			
            int total = 0;
			for (int i = 0; i < Groups.size(); i++) {
                int num = groupMap.getNumSeqs(Groups[i]);
                total += num;
				m->mothurOut(Groups[i] + " contains " + toString(num) + "."); m->mothurOutEndLine();
                out << Groups[i] << '\t' << num << endl;
			}
            out.close();
            m->mothurOut("\nTotal seqs: " + toString(total) + "."); m->mothurOutEndLine();
		}
        
        if (m->getControl_pressed()) { return 0; }
        
        if (countfile != "") {
            map<string, string> variables; 
            string thisOutputDir = outputDir;
            if (outputDir == "") {  thisOutputDir += util.hasPath(countfile);  }
            variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(countfile));
            string outputFileName = getOutputFileName("summary", variables);
            outputNames.push_back(outputFileName); outputTypes["summary"].push_back(outputFileName);
            ofstream out;
            util.openOutputFile(outputFileName, out);
            
			CountTable ct;
			ct.readTable(countfile, true, false);
        
			vector<string> nameGroups = ct.getNamesOfGroups();
			if (Groups.size() == 0) { Groups = nameGroups;  }
			
            int total = 0;
			for (int i = 0; i < Groups.size(); i++) {
                int num = ct.getGroupCount(Groups[i]);
                total += num;
				m->mothurOut(Groups[i] + " contains " + toString(num) + "."); m->mothurOutEndLine();
                out << Groups[i] << '\t' << num << endl;
			}
            out.close();
            
            m->mothurOut("\nTotal seqs: " + toString(total) + "."); m->mothurOutEndLine();
		}
		
		if (m->getControl_pressed()) { return 0; }
		
		if (sharedfile != "")		{		
			InputData input(sharedfile, "sharedfile", Groups);
			SharedRAbundVectors* lookup = input.getSharedRAbundVectors();
            Groups = lookup->getNamesGroups();
			
            map<string, string> variables; 
            string thisOutputDir = outputDir;
            if (outputDir == "") {  thisOutputDir += util.hasPath(sharedfile);  }
            variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(sharedfile));
            string outputFileName = getOutputFileName("summary", variables);
            outputNames.push_back(outputFileName); outputTypes["summary"].push_back(outputFileName);
            ofstream out;
            util.openOutputFile(outputFileName, out);
            
            int total = 0;
            vector<string> groups = lookup->getNamesGroups();
			for (int i = 0; i < groups.size(); i++) {
                int num = lookup->getNumSeqs(groups[i]);
                total += num;
				m->mothurOut(groups[i] + " contains " + toString(num) + "."); m->mothurOutEndLine();
                out << groups[i] << '\t' << num << endl;
			}
            out.close();
			
            m->mothurOut("\nTotal seqs: " + toString(total) + "."); m->mothurOutEndLine();
		}
			
        m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}	
		m->mothurOutEndLine();
        
		return 0;		
	}
	
	catch(exception& e) {
		m->errorOut(e, "CountGroupsCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************


