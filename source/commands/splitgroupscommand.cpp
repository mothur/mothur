/*
 *  splitgroupscommand.cpp
 *  Mothur
 *
 *  Created by westcott on 9/20/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "splitgroupscommand.h"

#include "sequenceparser.h"
#include "counttable.h"

//**********************************************************************************************************************
vector<string> SplitGroupCommand::setParameters(){	
	try {
        CommandParameter pflow("flow", "InputTypes", "", "", "none", "none", "none","fasta",false,false,true); parameters.push_back(pflow);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fasta",false,false,true); parameters.push_back(pfasta);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","name",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "CountGroup", "none","count",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "CountGroup", "none","group",false,false,true); parameters.push_back(pgroup);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitGroupCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SplitGroupCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The split.groups command reads a group or count file, and parses your files by groups. \n";
		helpString += "The split.groups command parameters are fasta, flow, name, group, count and groups.\n";
		helpString += "The group or count parameter is required.\n";
		helpString += "The groups parameter allows you to select groups to create files for.  \n";
		helpString += "For example if you set groups=A-B-C, you will get a .A.fasta, .A.names, .B.fasta, .B.names, .C.fasta, .C.names files.  \n";
		helpString += "If you want .fasta and .names files for all groups, set groups=all.  \n";
		helpString += "The split.groups command should be used in the following format: split.group(fasta=yourFasta, group=yourGroupFile).\n";
		helpString += "Example: split.groups(fasta=abrecovery.fasta, group=abrecovery.groups).\n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitGroupCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string SplitGroupCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta")        {  pattern = "[filename],[group],fasta";        }
        else if (type == "flow")    {  pattern = "[filename],[group],flow";         }
        else if (type == "name")    {  pattern = "[filename],[group],names";        }
        else if (type == "count")   {  pattern = "[filename],[group],count_table";  }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "SplitGroupCommand", "getOutputPattern");
        exit(1);
    }
}

//**********************************************************************************************************************
SplitGroupCommand::SplitGroupCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["name"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitGroupCommand", "SplitGroupCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
SplitGroupCommand::SplitGroupCommand(string option)  {
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
			outputTypes["fasta"] = tempOutNames;
            outputTypes["flow"] = tempOutNames;
			outputTypes["name"] = tempOutNames;
            outputTypes["count"] = tempOutNames;
		
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
				
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
                
                it = parameters.find("flow");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["flow"] = inputDir + it->second;		}
                }
			}

			
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") { namefile = ""; }	
			else { current->setNameFile(namefile); }
		
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not open") { abort = true; }
            else if (fastafile == "not found") {	fastafile = "";  }
            else { current->setFastaFile(fastafile); }
            
            flowfile = validParameter.validFile(parameters, "flow");
            if (flowfile == "not open") { abort = true; }
            else if (flowfile == "not found") {	flowfile = "";  }
            else { current->setFlowFile(flowfile); }
			
			groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not open") {  groupfile = ""; abort = true; }	
			else if (groupfile == "not found") { groupfile = "";
			}else {  current->setGroupFile(groupfile); }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { countfile = ""; abort = true; }
			else if (countfile == "not found") { countfile = ""; }	
			else { current->setCountFile(countfile); }
            
            if ((fastafile == "") && (flowfile == "")) {
                fastafile = current->getFastaFile();
                if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter.\n"); }
                else {
                    flowfile = current->getFlowFile();
                    if (flowfile != "") {  m->mothurOut("Using " + flowfile + " as input file for the flow parameter.\n");  }
                    else { m->mothurOut("[ERROR]: You need to provide a fasta or flow file.\n");  abort = true; }
                }
            }
            
            if ((countfile != "") && (namefile != "")) { m->mothurOut("You must enter ONLY ONE of the following: count or name."); m->mothurOutEndLine(); abort = true; }
            
            if ((countfile != "") && (groupfile != "")) { m->mothurOut("You must enter ONLY ONE of the following: count or group."); m->mothurOutEndLine(); abort = true; }
            
            if ((countfile == "") && (groupfile == "")) {
                if (namefile == "") { //check for count then group
                    countfile = current->getCountFile(); 
					if (countfile != "") {  m->mothurOut("Using " + countfile + " as input file for the count parameter."); m->mothurOutEndLine(); }
					else { 
						groupfile = current->getGroupFile(); 
                        if (groupfile != "") {  m->mothurOut("Using " + groupfile + " as input file for the group parameter."); m->mothurOutEndLine(); }
                        else { 
                            m->mothurOut("You need to provide a count or group file."); m->mothurOutEndLine(); 
                            abort = true; 
                        }	
					}	
                }else { //check for group
                    groupfile = current->getGroupFile(); 
                    if (groupfile != "") {  m->mothurOut("Using " + groupfile + " as input file for the group parameter."); m->mothurOutEndLine(); }
                    else { 
                        m->mothurOut("You need to provide a count or group file."); m->mothurOutEndLine(); 
                        abort = true; 
                    }	
                }
            }
			
			groups = validParameter.valid(parameters, "groups");		
			if (groups == "not found") { groups = ""; }
			else { util.splitAtDash(groups, Groups);
                    if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }	}
						
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	
                if (groupfile != "") { outputDir = util.hasPath(groupfile); }
                else { outputDir = util.hasPath(countfile);  }
            }
			
            if (countfile == "") {
                if (namefile == "") {
                    vector<string> files;
                    if (fastafile != "")    {  files.push_back(fastafile);  }
                    else                    {  files.push_back(flowfile);   }
                    if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
                }
            }
		}

	}
	catch(exception& e) {
		m->errorOut(e, "SplitGroupCommand", "SplitAbundCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
int SplitGroupCommand::execute(){
	try {
	
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
        if (flowfile != "")     {  splitFlow();  }
        if (fastafile != "")    {
            if (countfile == "" )   {  runNameGroup();  }
            else                    {  runCount();      }
        }
				
		if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);	} return 0; }
		
		string currentName = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); }
		}
        
        itTypes = outputTypes.find("flow");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFlowFile(currentName); }
        }
		
		itTypes = outputTypes.find("name");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setNameFile(currentName); }
		}
        
        itTypes = outputTypes.find("count");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setCountFile(currentName); }
		}
		
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitGroupCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int SplitGroupCommand::runNameGroup(){
	try {
        SequenceParser* parser;
		if (namefile == "") {	parser = new SequenceParser(groupfile, fastafile, Groups);				}
		else				{	parser = new SequenceParser(groupfile, fastafile, namefile, Groups);	}
		
		if (m->getControl_pressed()) { delete parser; return 0; }
        
		vector<string> namesGroups = parser->getNamesOfGroups();
        if (Groups.size() == 0) { Groups = namesGroups; }
		
		string fastafileRoot = outputDir + util.getRootName(util.getSimpleName(fastafile));
		string namefileRoot = outputDir + util.getRootName(util.getSimpleName(namefile));
		
		m->mothurOutEndLine();
		for (int i = 0; i < Groups.size(); i++) {
			
			m->mothurOut("Processing group: " + Groups[i] + "\n"); 
			
            map<string, string> variables; 
            variables["[filename]"] = fastafileRoot;
            variables["[group]"] = Groups[i];

			string newFasta = getOutputFileName("fasta",variables);
            variables["[filename]"] = namefileRoot;
			string newName = getOutputFileName("name",variables);
			
            long long numSeqs = 0;
			parser->getSeqs(Groups[i], newFasta, "/ab=", "/", numSeqs, false);
			outputNames.push_back(newFasta); outputTypes["fasta"].push_back(newFasta);
			if (m->getControl_pressed()) { delete parser; for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);	} return 0; }
            
			if (namefile != "") { 
				parser->getNameMap(Groups[i], newName); 
				outputNames.push_back(newName); outputTypes["name"].push_back(newName);
			}
			
			if (m->getControl_pressed()) { delete parser; for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);	} return 0; }
		}
		
		delete parser;
        
        return 0;

    }
	catch(exception& e) {
		m->errorOut(e, "SplitGroupCommand", "runNameGroup");
		exit(1);
	}
}
//**********************************************************************************************************************
struct flowOutput {
    string output;
    string filename;
    int total;
    
    flowOutput(string f) { filename = f; output = ""; total = 0;  }
    flowOutput() { filename = ""; output = ""; total = 0;  }
    flowOutput(string f, string o, int t) : filename(f), output(o), total(t) {}
    
};
//**********************************************************************************************************************
int SplitGroupCommand::splitFlow(){
    try {
        GroupMap* groupMap = NULL;
        CountTable* ct = NULL;
        vector<string> namesGroups;
        if (groupfile != "") {
            groupMap = new GroupMap(groupfile);
            groupMap->readMap();
            namesGroups = groupMap->getNamesOfGroups();
        }else if (countfile != ""){
            ct = new CountTable();
            ct->readTable(countfile, true, true);
            namesGroups = ct->getNamesOfGroups();
        }else { m->mothurOut("[ERROR]: you must provide a count or group file to split by group. quitting... \n"); m->setControl_pressed(true);  }
        
        if (Groups.size() == 0) { Groups = namesGroups; }
        
        if (m->getControl_pressed()) { if (groupMap != NULL) { delete groupMap; }else if (ct != NULL) { delete ct; } return 0; }
        
        string flowfileRoot = outputDir + util.getRootName(util.getSimpleName(flowfile));
        
        ifstream in; int numFlows = 0;
        util.openInputFile(flowfile, in);
        in >> numFlows; util.gobble(in);
        
        map<string, flowOutput> parsedFlowData;
        for (int i = 0; i < Groups.size(); i++) {
            map<string, string> variables;
            variables["[filename]"] = flowfileRoot;
            variables["[group]"] = Groups[i];
            string newFlow = getOutputFileName("flow",variables);
            
            flowOutput thisGroupsInfo(newFlow);
            parsedFlowData[Groups[i]] = thisGroupsInfo;
            
            ofstream out;
            util.openOutputFile(newFlow, out); out << numFlows << endl; out.close();  //clear file for append
        
            if (m->getControl_pressed()) { break; }
        }
        
        string name, flows;
        int count = 0;
        while (!in.eof()) {
            if (m->getControl_pressed()) { break; }
            
            in >> name; util.gobble(in);
            flows = util.getline(in); util.gobble(in);
            
            vector<string> thisSeqsGroups;
            if (groupMap != NULL) {
                string thisGroup = groupMap->getGroup(name);
                thisSeqsGroups.push_back(thisGroup);
            }else if (ct != NULL) { thisSeqsGroups  = ct->getGroups(name); }
            
            for (int i = 0; i < thisSeqsGroups.size(); i++) {
                
                map<string, flowOutput>::iterator it = parsedFlowData.find(thisSeqsGroups[i]);
                
                if (it != parsedFlowData.end()) {
                    it->second.total++; it->second.output += name + ' ' + flows + '\n';
                    if (it->second.total % 100 == 0) { //buffer write
                        ofstream out; util.openOutputFileAppend(it->second.filename, out);
                        out << it->second.output; it->second.output = "";
                    }
                } //else not in the groups we are looking to parse, so ignore
            }
            
            count++;
        }
        
        //output rest
        for (map<string, flowOutput>::iterator it = parsedFlowData.begin(); it != parsedFlowData.end(); it++) {
            if (m->getControl_pressed()) { break; }
            
            if (it->second.output != "") { //more seqs to output
                ofstream out; util.openOutputFileAppend(it->second.filename, out);
                out << it->second.output; it->second.output = ""; outputNames.push_back(it->second.filename); outputTypes["flow"].push_back(it->second.filename);
            }else if (it->second.total == 0) { //no seqs for this group, remove file
                util.mothurRemove(it->second.filename);
            }else { //finished writing, just add to list of output files
                outputNames.push_back(it->second.filename); outputTypes["flow"].push_back(it->second.filename);
            }
        }
        
        if (m->getControl_pressed()) { if (groupMap != NULL) { delete groupMap; }else if (ct != NULL) { delete ct; } return 0; }
        
        return count;
        
    }
    catch(exception& e) {
        m->errorOut(e, "SplitGroupCommand", "splitFlow");
        exit(1);
    }
}
//**********************************************************************************************************************
int SplitGroupCommand::runCount(){
	try {
        
        CountTable ct;
        ct.readTable(countfile, true, false);
        if (!ct.hasGroupInfo()) { m->mothurOut("[ERROR]: your count file does not contain group info, cannot split by group.\n"); m->setControl_pressed(true); }
        
        if (m->getControl_pressed()) { return 0; }
        
        vector<string> namesGroups = ct.getNamesOfGroups();
        
        //fill filehandles with neccessary ofstreams
        map<string, string> ffiles; //group -> filename
        map<string, string> cfiles; //group -> filename
        for (int i=0; i<Groups.size(); i++) {
            ofstream ftemp, ctemp;
            map<string, string> variables; 
            variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(fastafile));
            variables["[group]"] = Groups[i];
            string newFasta = getOutputFileName("fasta",variables);
            outputNames.push_back(newFasta); outputTypes["fasta"].push_back(newFasta);
            ffiles[Groups[i]] = newFasta;
            util.openOutputFile(newFasta, ftemp); ftemp.close();
            
            variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(countfile));
            string newCount = getOutputFileName("count",variables);
            outputNames.push_back(newCount); outputTypes["count"].push_back(newCount);
            cfiles[Groups[i]] = newCount;
            util.openOutputFile(newCount, ctemp);
            ctemp << "Representative_Sequence\ttotal\t" << Groups[i] << endl; ctemp.close();
        }
        
        ifstream in; 
        util.openInputFile(fastafile, in);
        
        while (!in.eof()) {
            Sequence seq(in); util.gobble(in);
            
            if (m->getControl_pressed()) { break; }
            if (seq.getName() != "") {
                vector<string> thisSeqsGroups = ct.getGroups(seq.getName());
                for (int i = 0; i < thisSeqsGroups.size(); i++) {
                    if (util.inUsersGroups(thisSeqsGroups[i], Groups)) { //if this sequence belongs to a group we want them print
                        ofstream outf, outc;
                        util.openOutputFileAppend(ffiles[thisSeqsGroups[i]], outf);
                        seq.printSequence(outf); outf.close();
                        int numSeqs = ct.getGroupCount(seq.getName(), thisSeqsGroups[i]);
                        util.openOutputFileAppend(cfiles[thisSeqsGroups[i]], outc);
                        outc << seq.getName() << '\t' << numSeqs << '\t' << numSeqs << endl; outc.close();
                    }
                }
            }
        }
        in.close();
        
        return 0;

    }
	catch(exception& e) {
		m->errorOut(e, "SplitGroupCommand", "runCount");
		exit(1);
	}
}
//**********************************************************************************************************************


