/*
 *  removegroupscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 11/10/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "removegroupscommand.h"
#include "sequence.hpp"
#include "listvector.hpp"

#include "inputdata.h"
#include "designmap.h"

//**********************************************************************************************************************
vector<string> RemoveGroupsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "FNGLT","fasta",false,false,true); parameters.push_back(pfasta);
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "sharedGroup", "none","shared",false,false,true); parameters.push_back(pshared);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","name",false,false,true); parameters.push_back(pname);
        CommandParameter pphylip("phylip", "InputTypes", "", "", "none", "PhylipColumn", "none","phylip",false,false,true); parameters.push_back(pphylip);
        CommandParameter pcolumn("column", "InputTypes", "", "", "none", "PhylipColumn", "none","column",false,false,true); parameters.push_back(pcolumn);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","count",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "sharedGroup", "FNGLT","group",false,false,true); parameters.push_back(pgroup);	        
        CommandParameter pdesign("design", "InputTypes", "", "", "none", "sharedGroup", "FNGLT","design",false,false); parameters.push_back(pdesign);
		CommandParameter plist("list", "InputTypes", "", "", "none", "none", "FNGLT","list",false,false,true); parameters.push_back(plist);
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "none", "FNGLT","taxonomy",false,false,true); parameters.push_back(ptaxonomy);
		CommandParameter paccnos("accnos", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(paccnos);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
        CommandParameter psets("sets", "String", "", "", "", "", "","",false,false); parameters.push_back(psets);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["fasta"] = tempOutNames;
        outputTypes["taxonomy"] = tempOutNames;
        outputTypes["name"] = tempOutNames;
        outputTypes["group"] = tempOutNames;
        outputTypes["list"] = tempOutNames;
        outputTypes["shared"] = tempOutNames;
        outputTypes["design"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
        outputTypes["phylip"] = tempOutNames;
        outputTypes["column"] = tempOutNames;
        
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string RemoveGroupsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The remove.groups command removes sequences from a specfic group or set of groups from the following file types: fasta, name, group, count, list, taxonomy, design, phylip, column or sharedfile.\n";
		helpString += "It outputs a file containing the sequences NOT in the those specified groups, or with a sharedfile eliminates the groups you selected.\n";
		helpString += "The remove.groups command parameters are accnos, fasta, name, group, list, taxonomy, shared, design, phylip, column, sets and groups. The group or count parameter is required, unless you have a current group or count file or are using a sharedfile.\n";
		helpString += "You must also provide an accnos containing the list of groups to remove or set the groups or sets parameter to the groups you wish to remove.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like removed.  You can separate group names with dashes.\n";
        helpString += "The sets parameter allows you to specify which of the sets in your designfile you would like to remove.  You can separate set names with dashes.\n";
		helpString += "The remove.groups command should be in the following format: remove.groups(accnos=yourAccnos, fasta=yourFasta, group=yourGroupFile).\n";
		helpString += "Example remove.groups(accnos=amazon.accnos, fasta=amazon.fasta, group=amazon.groups).\n";
		helpString += "or remove.groups(groups=pasture, fasta=amazon.fasta, amazon.groups).\n";
		;
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string RemoveGroupsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta")            {   pattern = "[filename],pick,[extension]";    }
        else if (type == "taxonomy")    {   pattern = "[filename],pick,[extension]";    }
        else if (type == "name")        {   pattern = "[filename],pick,[extension]";    }
        else if (type == "group")       {   pattern = "[filename],pick,[extension]";    }
        else if (type == "count")       {   pattern = "[filename],pick,[extension]";    }
        else if (type == "phylip")      {   pattern = "[filename],pick,[extension]";    }
        else if (type == "column")      {   pattern = "[filename],pick,[extension]";    }
        else if (type == "list")        {   pattern = "[filename],[tag],pick,[extension]";    }
        else if (type == "shared")      {   pattern = "[filename],[tag],pick,[extension]";    }
        else if (type == "design")      {   pattern = "[filename],[tag],pick,[extension]-[filename],pick,[extension]";    }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "RemoveGroupsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
RemoveGroupsCommand::RemoveGroupsCommand(string option) : Command()  {
	try {

		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
            
			
			//check for required parameters
			accnosfile = validParameter.validFile(parameters, "accnos");
			if (accnosfile == "not open") { accnosfile = ""; abort = true; }
			else if (accnosfile == "not found") {  accnosfile = ""; }	
			else { current->setAccnosFile(accnosfile); }
            
            phylipfile = validParameter.validFile(parameters, "phylip");
            if (phylipfile == "not open") { phylipfile = ""; abort = true; }
            else if (phylipfile == "not found") { phylipfile = ""; }
            else { 	current->setPhylipFile(phylipfile); }
            
            columnfile = validParameter.validFile(parameters, "column");
            if (columnfile == "not open") { columnfile = ""; abort = true; }
            else if (columnfile == "not found") { columnfile = ""; }
            else {  current->setColumnFile(columnfile);	}
			
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not open") { fastafile = ""; abort = true; }
			else if (fastafile == "not found") {  fastafile = "";  }	
			else { current->setFastaFile(fastafile); }
			
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") {  namefile = "";  }	
			else { current->setNameFile(namefile); }
			
			groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not open") { groupfile = "";  abort = true; }
			else if (groupfile == "not found") {  	groupfile = "";		}
			else { current->setGroupFile(groupfile); }	
			
			listfile = validParameter.validFile(parameters, "list");
			if (listfile == "not open") { listfile = ""; abort = true; }
			else if (listfile == "not found") {  listfile = "";  }
			else { current->setListFile(listfile); }
			
			taxfile = validParameter.validFile(parameters, "taxonomy");
			if (taxfile == "not open") { taxfile = ""; abort = true; }
			else if (taxfile == "not found") {  taxfile = "";  }
			else { current->setTaxonomyFile(taxfile); }
            
            designfile = validParameter.validFile(parameters, "design");
			if (designfile == "not open") { designfile = ""; abort = true; }
			else if (designfile == "not found") {  designfile = "";  }
			else { current->setDesignFile(designfile); }
			
			groups = validParameter.valid(parameters, "groups");
            if (groups == "not found") { groups = ""; }
            else { util.splitAtDash(groups, Groups);  }
            
            sets = validParameter.valid(parameters, "sets");
            if (sets == "not found") { sets = ""; }
            else { util.splitAtDash(sets, Sets);  }
			
			sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }
			else if (sharedfile == "not found") {  sharedfile = "";  }
			else { current->setSharedFile(sharedfile); }
			
			
			countfile = validParameter.validFile(parameters, "count");
            if (countfile == "not open") { countfile = ""; abort = true; }
            else if (countfile == "not found") { countfile = "";  }	
            else { current->setCountFile(countfile); }
            
            if ((namefile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: name or count.\n");  abort = true;
            }
            
            if ((groupfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: group or count.\n");  abort=true;
            }
            
			
			if ((sharedfile == "") && (groupfile == "") && (designfile == "") && (countfile == "")) { 
				//is there are current file available for any of these?
				if ((namefile != "") || (fastafile != "") || (listfile != "") || (taxfile != "")) {
					//give priority to group, then shared
					groupfile = current->getGroupFile(); 
					if (groupfile != "") {  m->mothurOut("Using " + groupfile + " as input file for the group parameter.\n");  }
					else { 
						sharedfile = current->getSharedFile(); 
						if (sharedfile != "") { m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.\n");  }
						else { 
							countfile = current->getCountFile(); 
                            if (countfile != "") { m->mothurOut("Using " + countfile + " as input file for the count parameter.\n");  }
                            else { 
                                m->mothurOut("[ERROR]: You have no current groupfile, countfile or sharedfile and one is required.\n");  abort = true;
                            }
						}
					}
				}else {
					//give priority to shared, then group
					sharedfile = current->getSharedFile(); 
					if (sharedfile != "") {  m->mothurOut("Using " + sharedfile + " as input file for the shared parameter.\n");  }
					else { 
						groupfile = current->getGroupFile(); 
						if (groupfile != "") { m->mothurOut("Using " + groupfile + " as input file for the group parameter.\n");  }
						else { 
							designfile = current->getDesignFile(); 
                            if (designfile != "") { m->mothurOut("Using " + designfile + " as input file for the design parameter.\n");  }
                            else { 
                                countfile = current->getCountFile(); 
                                if (countfile != "") { m->mothurOut("Using " + countfile + " as input file for the count parameter.\n");  }
                                else { 
                                    m->mothurOut("[ERROR]: You have no current groupfile, designfile, countfile or sharedfile and one is required.\n");  abort = true;
                                }
                                
                            }
						}
					}
				}
			}
			
			if ((accnosfile == "") && (Groups.size() == 0) && (Sets.size() == 0)) { m->mothurOut("[ERROR]: You must provide an accnos file or specify groups using the groups or sets parameters.\n"); abort = true; }
            
            if ((Groups.size() != 0) && (Sets.size() != 0)) { m->mothurOut("[ERROR]: You cannot use the groups and sets parameters at the same time, quitting.\n"); abort = true; }
            
            if ((Sets.size() != 0) && (designfile == "")) { m->mothurOut("[ERROR]: You must provide a design file when using the sets parameter.\n"); abort = true;  }
            
			if ((phylipfile == "") && (columnfile == "") && (fastafile == "") && (namefile == "") && (countfile == "") && (groupfile == "")  && (designfile == "") && (sharedfile == "") && (listfile == "") && (taxfile == ""))  { m->mothurOut("[ERROR]: You must provide at least one of the following: fasta, name, taxonomy, group, shared, design, count, phylip, column or list.\n");  abort = true; }
			if (((groupfile == "") && (countfile == "")) && ((namefile != "") || (fastafile != "") || (listfile != "") || (taxfile != "")))  { m->mothurOut("[ERROR]: If using a fasta, name, taxonomy, group or list, then you must provide a group or count file.\n");  abort = true; }
            
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "RemoveGroupsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int RemoveGroupsCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		//get groups you want to remove
		if (accnosfile != "") { util.readAccnos(accnosfile, Groups);   }
        else if (Sets.size() != 0) {  fillGroupsFromDesign(); }
		
		if (groupfile != "") {
			groupMap = new GroupMap(groupfile);
			groupMap->readMap();
			
			//make sure groups are valid
			//takes care of user setting groupNames that are invalid or setting groups=all
			vector<string> namesGroups = groupMap->getNamesOfGroups();
			vector<string> checkedGroups;
            for (int i = 0; i < Groups.size(); i++) {
                if (util.inUsersGroups(Groups[i], namesGroups)) { checkedGroups.push_back(Groups[i]); }
                else {  m->mothurOut("[WARNING]: " + Groups[i] + " is not a valid group in your groupfile, ignoring.\n"); }
            }
            
            if (checkedGroups.size() == 0) { m->mothurOut("[ERROR]: no valid groups, aborting.\n"); delete groupMap; return 0; }
			else { Groups = checkedGroups; }
            
			//fill names with names of sequences that are from the groups we want to remove 
			fillNames();
			
			delete groupMap;
		}else if (countfile != ""){
            if ((fastafile != "") || (listfile != "") || (taxfile != "")) { 
                //m->mothurOut("\n[NOTE]: The count file should contain only unique names, so mothur assumes your fasta, list and taxonomy files also contain only uniques.\n\n");
            }
            vector<string> gNamesOfGroups;
            CountTable ct;
            
            if (!ct.testGroups(countfile, gNamesOfGroups)) { m->mothurOut("[ERROR]: your count file does not contain group info, aborting.\n"); return 0; }
            
            if (Groups.size() == 0) { return 0; }
            
            ct.readTable(countfile, true, false);
            int oldTotal = ct.getNumSeqs();
            
            vector<string> namesOfSeqs = ct.getNamesOfSeqs(); //all names
            
            for (int i = 0; i < Groups.size(); i++) { ct.removeGroup(Groups[i]); }
            
            vector<string> newSeqs = ct.getNamesOfSeqs(); //names of seqs left after removing groups
            
            set<string> goodNames = util.mothurConvert(newSeqs);
            
            for (int i = 0; i < namesOfSeqs.size(); i++) {
                if (goodNames.count(namesOfSeqs[i]) == 0) { //you aren't on good list
                    names.insert(namesOfSeqs[i]); //add to remove list
                }
            }
            
            string thisOutputDir = outputdir;
            if (outputdir == "") {  thisOutputDir += util.hasPath(countfile);  }
            map<string, string> variables;
            variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(countfile));
            variables["[extension]"] = util.getExtension(countfile);
            string outputFileName = getOutputFileName("count", variables);
            
            int selectedCount = ct.getNumSeqs();
            
            if (selectedCount == 0) {  m->mothurOut("Your file does NOT contain sequences from the groups you wish to get.\n");   }
            else {
                ct.printTable(outputFileName);
                outputTypes["count"].push_back(outputFileName); outputNames.push_back(outputFileName);
            }
            
            int removedCount = oldTotal - selectedCount;
            
            m->mothurOut("Removed " + toString(removedCount) + " sequences from your count file.\n");
        }

				
		if (m->getControl_pressed()) { return 0; }
		
		//read through the correct file and output lines you want to keep
		if (namefile != "")			{		readName();		}
		if (fastafile != "")		{		readFasta();	}
		if (groupfile != "")		{		readGroup();	}
		if (listfile != "")			{		readList();		}
		if (taxfile != "")			{		readTax();		}
		if (sharedfile != "")		{		readShared();	}
        if (designfile != "")		{		readDesign();	}
        if (phylipfile != "")		{		readPhylip();		}
        if (columnfile != "")		{		readColumn();       }
		
		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); } return 0; }
				
		if (outputNames.size() != 0) {
			m->mothurOutEndLine();
			m->mothurOut("Output File names: \n");
			for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
			m->mothurOutEndLine();
			
			//set fasta file as new current fastafile
			string currentName = "";
			itTypes = outputTypes.find("fasta");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); }
			}
			
			itTypes = outputTypes.find("name");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setNameFile(currentName); }
			}
			
			itTypes = outputTypes.find("group");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setGroupFile(currentName); }
			}
			
			itTypes = outputTypes.find("list");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setListFile(currentName); }
			}
			
			itTypes = outputTypes.find("taxonomy");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setTaxonomyFile(currentName); }
			}
			
			itTypes = outputTypes.find("shared");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setSharedFile(currentName); }
			}
            
            itTypes = outputTypes.find("design");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setDesignFile(currentName); }
			}
            
            itTypes = outputTypes.find("count");
			if (itTypes != outputTypes.end()) {
				if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setCountFile(currentName); }
			}
            
            itTypes = outputTypes.find("phylip");
            if (itTypes != outputTypes.end()) {
                if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setPhylipFile(currentName); }
            }
            
            itTypes = outputTypes.find("column");
            if (itTypes != outputTypes.end()) {
                if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setColumnFile(currentName); }
            }
		}
		
		return 0;		
	}
	
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "execute");
		exit(1);
	}
}

//**********************************************************************************************************************
void RemoveGroupsCommand::readFasta(){
	try {
		string thisOutputDir = outputdir;
		if (outputdir == "") {  thisOutputDir += util.hasPath(fastafile);  }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(fastafile));
        variables["[extension]"] = util.getExtension(fastafile);
		string outputFileName = getOutputFileName("fasta", variables);
		
		ofstream out; util.openOutputFile(outputFileName, out);
		ifstream in; util.openInputFile(fastafile, in);
		string name;
		
		bool wroteSomething = false;
		int removedCount = 0;
		
		while(!in.eof()){
			if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return; }
			
			Sequence currSeq(in);
			name = currSeq.getName();
			
			if (name != "") {
				//if this name is in the accnos file
				if (names.count(name) == 0) {
					wroteSomething = true;
					currSeq.printSequence(out); 
				}else { 
					//if you are not in the accnos file check if you are a name that needs to be changed
					map<string, string>::iterator it = uniqueToRedundant.find(name);
					if (it != uniqueToRedundant.end()) {
						wroteSomething = true;
						currSeq.setName(it->second);
						currSeq.printSequence(out);
					}else { removedCount++; }
				}
			}
			util.gobble(in);
		}
		in.close();	
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the groups you wish to remove.\n");   }
		outputTypes["fasta"].push_back(outputFileName);  outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from your fasta file.\n");
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "readFasta");
		exit(1);
	}
}
//**********************************************************************************************************************
void RemoveGroupsCommand::readShared(){
	try {
		string thisOutputDir = outputdir;
		if (outputdir == "") {  thisOutputDir += util.hasPath(sharedfile);  }
		
		//get group names from sharedfile so we can set Groups to the groupNames we want to keep
		//that way we can take advantage of the reads in inputdata and sharedRabundVector
		InputData input(sharedfile, "sharedfile", nullVector);
		SharedRAbundVectors* lookup = input.getSharedRAbundVectors();
        vector<string> allGroupsNames = lookup->getNamesGroups();
        
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(sharedfile));
        variables["[extension]"] = util.getExtension(sharedfile);
	
		vector<string> groupsToKeep;
		for (int i = 0; i < allGroupsNames.size(); i++) {
			if (!util.inUsersGroups(allGroupsNames[i], Groups)) {
				groupsToKeep.push_back(allGroupsNames[i]);
			}
		}
		
        if (allGroupsNames.size() == groupsToKeep.size()) { m->mothurOut("Your file does not contain any groups you wish to remove.\n");  delete lookup;  return; }
		
		//reset read 
		delete lookup;
		
		InputData input2(sharedfile, "sharedfile", groupsToKeep);
		lookup = input2.getSharedRAbundVectors();

		bool wroteSomething = false;
        bool printHeaders = true;
		
		while(lookup != nullptr) {
			
			variables["[tag]"] = lookup->getLabel();
            string outputFileName = getOutputFileName("shared", variables);
			ofstream out;
			util.openOutputFile(outputFileName, out);
			outputTypes["shared"].push_back(outputFileName);  outputNames.push_back(outputFileName);
			
            if (m->getControl_pressed()) { out.close();  util.mothurRemove(outputFileName);  delete lookup; return; }
			
            if (lookup->size() != 0) { wroteSomething = true; }
            lookup->print(out, printHeaders);
			
			//get next line to process
			//prevent memory leak
            delete lookup;
			lookup = input2.getSharedRAbundVectors();
			
			out.close();
		}
        
		if (wroteSomething == false) {  m->mothurOut("Your file contains only the groups you wish to remove.\n");   }
		
		string groupsString = "";
		for (int i = 0; i < Groups.size()-1; i++) {	groupsString += Groups[i] + ", "; }
		groupsString += Groups[Groups.size()-1];
		
		m->mothurOut("Removed groups: " + groupsString + " from your shared file.\n");
    }
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "readShared");
		exit(1);
	}
}
//**********************************************************************************************************************
void RemoveGroupsCommand::readList(){
	try {
		string thisOutputDir = outputdir;
		if (outputdir == "") {  thisOutputDir += util.hasPath(listfile);  }
		map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(listfile));
        variables["[extension]"] = util.getExtension(listfile);
        InputData input(listfile, "list", nullVector);
        ListVector* list = input.getListVector();
				
        bool wroteSomething = false;
		int removedCount = 0;
		
		while(list != nullptr) {
			
			removedCount = 0;
			
            variables["[tag]"] = list->getLabel();
            string outputFileName = getOutputFileName("list", variables);
			
			ofstream out;
			util.openOutputFile(outputFileName, out);
			outputTypes["list"].push_back(outputFileName);  outputNames.push_back(outputFileName);
            
            vector<string> binLabels = list->getLabels();
            vector<string> newBinLabels;
			
			//make a new list vector
			ListVector newList;
			newList.setLabel(list->getLabel());
			
			//for each bin
			for (int i = 0; i < list->getNumBins(); i++) {
				if (m->getControl_pressed()) {  out.close();  util.mothurRemove(outputFileName);  return; }
				
				//parse out names that are in accnos file
				string binnames = list->get(i);
				
				string newNames = "";
				while (binnames.find_first_of(',') != -1) { 
					string name = binnames.substr(0,binnames.find_first_of(','));
					binnames = binnames.substr(binnames.find_first_of(',')+1, binnames.length());
					
					//if that name is in the .accnos file, add it
					if (names.count(name) == 0) {  newNames += name + ",";  }
					else {
						//if you are not in the accnos file check if you are a name that needs to be changed
						map<string, string>::iterator it = uniqueToRedundant.find(name);
						if (it != uniqueToRedundant.end()) {
							newNames += it->second + ",";
						}else { removedCount++; }
					}
				}
				
				//get last name
				if (names.count(binnames) == 0) {  newNames += binnames + ",";  }
				else { //if you are not in the accnos file check if you are a name that needs to be changed
					map<string, string>::iterator it = uniqueToRedundant.find(binnames);
					if (it != uniqueToRedundant.end()) {
						newNames += it->second + ",";
					}else { removedCount++; }
				}
				
				//if there are names in this bin add to new list
				if (newNames != "") {  
					newNames = newNames.substr(0, newNames.length()-1); //rip off extra comma
					newList.push_back(newNames);
                    newBinLabels.push_back(binLabels[i]);
				}
			}
			
			//print new listvector
			if (newList.getNumBins() != 0) {
				wroteSomething = true;
				newList.setLabels(newBinLabels);
                newList.print(out, false);
			}
			
            out.close();
            
            delete list;
            list = input.getListVector();
		}
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the groups you wish to remove.\n");   }
		m->mothurOut("Removed " + toString(removedCount) + " sequences from your list file.\n");
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "readList");
		exit(1);
	}
}
//**********************************************************************************************************************
void RemoveGroupsCommand::readName(){
	try {
		string thisOutputDir = outputdir;
		if (outputdir == "") {  thisOutputDir += util.hasPath(namefile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(namefile));
        variables["[extension]"] = util.getExtension(namefile);
		string outputFileName = getOutputFileName("name", variables);	
		
        ofstream out; util.openOutputFile(outputFileName, out);
		ifstream in; util.openInputFile(namefile, in);
		string name, firstCol, secondCol;
		
		bool wroteSomething = false;
		int removedCount = 0;
		
		while(!in.eof()){
			if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return; }
			
			in >> firstCol;		util.gobble(in);		
			in >> secondCol;	util.gobble(in);
			
			vector<string> parsedNames;
			util.splitAtComma(secondCol, parsedNames);
						
			vector<string> validSecond;  validSecond.clear();
			for (int i = 0; i < parsedNames.size(); i++) {
				if (names.count(parsedNames[i]) == 0) {
					validSecond.push_back(parsedNames[i]);
				}
			}
			
			removedCount += parsedNames.size()-validSecond.size();
			
			//if the name in the first column is in the set then print it and any other names in second column also in set
			if (names.count(firstCol) == 0) {
				
				wroteSomething = true;
				
				out << firstCol << '\t';
				
				//you know you have at least one valid second since first column is valid
				for (int i = 0; i < validSecond.size()-1; i++) {  out << validSecond[i] << ',';  }
				out << validSecond[validSecond.size()-1] << endl;
				
				//make first name in set you come to first column and then add the remaining names to second column
			}else {
				
				//you want part of this row
				if (validSecond.size() != 0) {
					
					wroteSomething = true;
					
					out << validSecond[0] << '\t';
					
					//you know you have at least one valid second since first column is valid
					for (int i = 0; i < validSecond.size()-1; i++) {  out << validSecond[i] << ',';  }
					out << validSecond[validSecond.size()-1] << endl;
					uniqueToRedundant[firstCol] = validSecond[0];
				}
			}
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the groups you wish to remove.\n");   }
		outputTypes["name"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from your name file.\n");
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "readName");
		exit(1);
	}
}
//**********************************************************************************************************************
void RemoveGroupsCommand::readGroup(){
	try {
		string thisOutputDir = outputdir;
		if (outputdir == "") {  thisOutputDir += util.hasPath(groupfile);  }
        map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(groupfile));
        variables["[extension]"] = util.getExtension(groupfile);
		string outputFileName = getOutputFileName("group", variables);	
		
        ofstream out; util.openOutputFile(outputFileName, out);
		ifstream in; util.openInputFile(groupfile, in);
		string name, group;
		
		bool wroteSomething = false;
		int removedCount = 0;
		
		while(!in.eof()){
			if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return; }
			
			in >> name;				//read from first column
			in >> group;			//read from second column
			
			//if this name is in the accnos file
			if (names.count(name) == 0) {
				wroteSomething = true;
				out << name << '\t' << group << endl;
			}else {  removedCount++;  }
			
			util.gobble(in);
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the groups you wish to remove.\n");   }
		outputTypes["group"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from your group file.\n");
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "readGroup");
		exit(1);
	}
}
//**********************************************************************************************************************
void RemoveGroupsCommand::readDesign(){
	try {
		string thisOutputDir = outputdir;
		if (outputdir == "") {  thisOutputDir += util.hasPath(designfile);  }
        map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(designfile));
        variables["[extension]"] = util.getExtension(designfile);
		string outputFileName = getOutputFileName("design", variables);
		
        DesignMap designMap(designfile);  if (m->getControl_pressed()) {  return ; }
        
        vector<string> groupsToKeep;
        vector<string> allGroups = designMap.getNamesGroups();
        for(int i = 0; i < allGroups.size(); i++) {
            if (!util.inUsersGroups(allGroups[i], Groups)) {
                groupsToKeep.push_back(allGroups[i]);
            }
        }
        
        bool wroteSomething = false;
        
        ofstream out;
        util.openOutputFile(outputFileName, out);
        
        int numGroupsFound = designMap.printGroups(out, groupsToKeep);
        
        if (numGroupsFound > 0) { wroteSomething = true; }
        
        out.close();
		
        names.clear(); names = util.mothurConvert(Groups);
        
        int removedCount = allGroups.size() - numGroupsFound;
        
		if (wroteSomething == false) {  m->mothurOut("Your file contains only groups from the groups you wish to remove.\n");   }
		outputTypes["design"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " groups from your design file.\n");
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "readDesign");
		exit(1);
	}
}

//**********************************************************************************************************************
void RemoveGroupsCommand::readTax(){
	try {
		string thisOutputDir = outputdir;
		if (outputdir == "") {  thisOutputDir += util.hasPath(taxfile);  }
		map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(taxfile));
        variables["[extension]"] = util.getExtension(taxfile);
		string outputFileName = getOutputFileName("taxonomy", variables);
		
        ofstream out; util.openOutputFile(outputFileName, out);
		ifstream in; util.openInputFile(taxfile, in);
		string name, tax;
		
		bool wroteSomething = false;
		int removedCount = 0;
		
		while(!in.eof()){
			if (m->getControl_pressed()) { in.close();  out.close();  util.mothurRemove(outputFileName);  return; }
			
            in >> name; util.gobble(in);
            tax = util.getline(in); util.gobble(in);
			
			//if this name is in the accnos file
			if (names.count(name) == 0) {
				wroteSomething = true;
				out << name << '\t' << tax << endl;
			}else {  //if you are not in the accnos file check if you are a name that needs to be changed
				map<string, string>::iterator it = uniqueToRedundant.find(name);
				if (it != uniqueToRedundant.end()) {
					wroteSomething = true;
					out << it->second << '\t' << tax << endl;
				}else { removedCount++; }  }
		}
		in.close();
		out.close();
		
		if (wroteSomething == false) {  m->mothurOut("Your file contains only sequences from the groups you wish to remove.\n");   }
		outputTypes["taxonomy"].push_back(outputFileName); outputNames.push_back(outputFileName);
		
		m->mothurOut("Removed " + toString(removedCount) + " sequences from your taxonomy file.\n");
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "readTax");
		exit(1);
	}
}
//**********************************************************************************************************************
void RemoveGroupsCommand::readPhylip(){
    try {
        string thisOutputDir = outputdir;
        if (outputdir == "") {  thisOutputDir += util.hasPath(phylipfile);  }
        map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(phylipfile));
        variables["[extension]"] = util.getExtension(phylipfile);
        string outputFileName = getOutputFileName("phylip", variables);
        
        ifstream in; util.openInputFile(phylipfile, in);
        
        float distance;
        int square, nseqs; square = 0;
        string name;
        unsigned int row;
        set<unsigned int> rows; //converts names in names to a index
        row = 0;
        
        string numTest;
        in >> numTest >> name;
        
        if (!util.isContainingOnlyDigits(numTest)) { m->mothurOut("[ERROR]: expected a number and got " + numTest + ", quitting.\n");  m->setControl_pressed(true); return; }
        else { convert(numTest, nseqs); }
        
        //not one we want to remove
        if (names.count(name) == 0) { rows.insert(row); }
        row++;
        
        //is the matrix square?
        char d;
        while((d=in.get()) != EOF){
            
            if(isalnum(d)){
                square = 1;
                in.putback(d);
                for(int i=0;i<nseqs;i++){
                    in >> distance;
                }
                break;
            }
            if(d == '\n'){
                square = 0;
                break;
            }
        }
        
        //map name to row/column
        if(square == 0){
            for(int i=1;i<nseqs;i++){
                in >> name;
                if (names.count(name) == 0) { rows.insert(row); }
                row++;
                
                for(int j=0;j<i;j++){
                    if (m->getControl_pressed()) {  in.close(); return;  }
                    in >> distance;
                }
            }
        }
        else{
            for(int i=1;i<nseqs;i++){
                in >> name;
                if (names.count(name) == 0) { rows.insert(row);  }
                row++;
                for(int j=0;j<nseqs;j++){
                    if (m->getControl_pressed()) {  in.close(); return;  }
                    in >> distance;
                }
            }
        }
        in.close();
        
        if (m->getControl_pressed()) {  return; }
        
        //read through file only printing rows and columns of seqs in names
        ifstream inPhylip;
        util.openInputFile(phylipfile, inPhylip);
        
        inPhylip >> numTest;
        
        ofstream out;
        util.openOutputFile(outputFileName, out);
        outputTypes["phylip"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        out << (nseqs-names.size()) << endl;
        
        unsigned int count = 0;
        unsigned int keptCount = 0;
        if(square == 0){
            for(int i=0;i<nseqs;i++){
                inPhylip >> name;
                bool ignoreRow = false;
                
                if (names.count(name) != 0) { ignoreRow = true; count++; }
                else{ out << name; keptCount++; }
                
                for(int j=0;j<i;j++){
                    if (m->getControl_pressed()) {  inPhylip.close(); out.close();  return;  }
                    inPhylip >> distance;
                    if (!ignoreRow) {
                        //is this a column we want
                        if(rows.count(j) != 0) {  out << '\t' << distance;  }
                    }
                }
                if (!ignoreRow) { out << endl; }
            }
        }
        else{
            for(int i=0;i<nseqs;i++){
                inPhylip >> name;
                
                bool ignoreRow = false;
                
                if (names.count(name) != 0) { ignoreRow = true; count++; }
                else{ out << name; keptCount++; }
                
                for(int j=0;j<nseqs;j++){
                    if (m->getControl_pressed()) {  inPhylip.close(); out.close(); return;  }
                    inPhylip >> distance;
                    if (!ignoreRow) {
                        //is this a column we want
                        if(rows.count(j) != 0) {  out << '\t' << distance;  }
                    }
                }
                if (!ignoreRow) { out << endl; }
            }
        }
        inPhylip.close();
        out.close();
        
        if (keptCount == 0) {  m->mothurOut("Your file contains ONLY distances related to groups or sequences listed in the accnos file.\n");   }
        else if (count != names.size()) {
            m->mothurOut("[WARNING]: Your accnos file contains " + toString(names.size()) + " groups or sequences, but I only found " + toString(count) + " of them in the phylip file.\n");
            //rewrite with new number
            util.renameFile(outputFileName, outputFileName+".temp");
            ofstream out2;
            util.openOutputFile(outputFileName, out2);
            out2 << keptCount << endl;
            
            ifstream in3;
            util.openInputFile(outputFileName+".temp", in3);
            in3 >> nseqs; util.gobble(in3);
            char buffer[4096];
            while (!in3.eof()) {
                in3.read(buffer, 4096);
                out2.write(buffer, in3.gcount());
            }
            in3.close();
            out2.close();
            util.mothurRemove(outputFileName+".temp");
        }
        
        m->mothurOut("Removed " + toString(count) + " groups or sequences from your phylip file.\n");
    }
    catch(exception& e) {
        m->errorOut(e, "RemoveGroupsCommand", "readPhylip");
        exit(1);
    }
}
//**********************************************************************************************************************
void RemoveGroupsCommand::readColumn(){
    try {
        string thisOutputDir = outputdir;
        if (outputdir == "") {  thisOutputDir += util.hasPath(columnfile);  }
        map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(columnfile));
        variables["[extension]"] = util.getExtension(columnfile);
        string outputFileName = getOutputFileName("column", variables);
        outputTypes["column"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        
        ofstream out; util.openOutputFile(outputFileName, out);
        ifstream in; util.openInputFile(columnfile, in);
        
        set<string> removeNames;
        string firstName, secondName;
        float distance;
        bool wrote = false;
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { out.close(); in.close(); return; }
            
            in >> firstName >> secondName >> distance; util.gobble(in);
            
            //is either names in the accnos file
            if (names.count(firstName) != 0)       {
                removeNames.insert(firstName);
                if (names.count(secondName) != 0)  { removeNames.insert(secondName);      }   }
            else if (names.count(secondName) != 0) {
                removeNames.insert(secondName);
                if (names.count(firstName) != 0)   { removeNames.insert(firstName);     }   }
            else {
                wrote = true;
                out << firstName << '\t' << secondName << '\t' << distance << endl;
            }
        }
        in.close();
        out.close();
        
        if (!wrote) {  m->mothurOut("Your file contains ONLY distances related to groups or sequences listed in the accnos file.\n");   }
        else if (removeNames.size() != names.size()) {
            m->mothurOut("[WARNING]: Your accnos file contains " + toString(names.size()) + " groups or sequences, but I only found " + toString(removeNames.size()) + " of them in the column file.\n");
        }
        
        m->mothurOut("Removed " + toString(removeNames.size()) + " groups or sequences from your column file.\n");
    }
    catch(exception& e) {
        m->errorOut(e, "RemoveGroupsCommand", "readColumn");
        exit(1);
    }
}
//**********************************************************************************************************************
//names to remove
void RemoveGroupsCommand::fillNames(){
	try {
		vector<string> seqs = groupMap->getNamesSeqs();
		
		for (int i = 0; i < seqs.size(); i++) {
			
            if (m->getControl_pressed()) { break; }
			
			string group = groupMap->getGroup(seqs[i]);
			
			if (util.inUsersGroups(group, Groups)) {
				names.insert(seqs[i]);
			}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "RemoveGroupsCommand", "fillNames");
		exit(1);
	}
}
//**********************************************************************************************************************
void RemoveGroupsCommand::fillGroupsFromDesign(){
    try {
        DesignMap designMap(designfile);  if (m->getControl_pressed()) {  return ; }
        Groups = designMap.getNamesGroups(Sets);
    }
    catch(exception& e) {
        m->errorOut(e, "RemoveGroupsCommand", "fillGroupsFromDesign");
        exit(1);
    }
}
//**********************************************************************************************************************



