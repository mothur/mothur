/*
 *  splitabundcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 5/17/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "splitabundcommand.h"
#include "getseqscommand.h"
#include "getotuscommand.h"

//**********************************************************************************************************************
vector<string> SplitAbundCommand::setParameters(){	
	try {		
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fasta",false,false,true); parameters.push_back(pfasta);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "FNGLT", "none","name",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","count",false,false); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","group",false,false); parameters.push_back(pgroup);
		CommandParameter plist("list", "InputTypes", "", "", "none", "FNGLT", "none","list",false,false,true); parameters.push_back(plist);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pcutoff("cutoff", "Number", "", "0", "", "", "","",false,true); parameters.push_back(pcutoff);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        vector<string> tempOutNames;
        outputTypes["list"] = tempOutNames;
        outputTypes["name"] = tempOutNames;
        outputTypes["accnos"] = tempOutNames;
        outputTypes["group"] = tempOutNames;
        outputTypes["fasta"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
        
        abort = false; calledHelp = false;    allLines = true;
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SplitAbundCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The split.abund command reads a fasta file or a list or a names or a count file splits the sequences into rare and abundant groups. \n";
		helpString += "The split.abund command parameters are fasta, list, name, count, cutoff, group, label and cutoff.\n";
		helpString += "The fasta or a list or name or count parameter are required, and you must provide a cutoff value.\n";
		helpString += "The cutoff parameter is used to qualify what is abundant and rare. If cutoff < 1, mothur assumes this is a percentage. 0.02 -> rare reads represent <= 2% of total reads. \n";
		helpString += "The group parameter allows you to parse a group file into rare and abundant groups.\n";
		helpString += "The label parameter is used to read specific labels in your listfile you want to use.\n";
		helpString += "For example if you set groups=A-B-C, you will get a .A.abund, .A.rare, .B.abund, .B.rare, .C.abund, .C.rare files.  \n";
		helpString += "If you want .abund and .rare files for all groups, set groups=all.  \n";
		helpString += "The split.abund command should be used in the following format: split.abund(fasta=yourFasta, list=yourListFile, group=yourGroupFile, label=yourLabels, cutoff=yourCutoff).\n";
		helpString += "Example: split.abund(fasta=abrecovery.fasta, list=abrecovery.fn.list, group=abrecovery.groups, label=0.03, cutoff=2).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "getHelpString");
		exit(1);
	}
}

//**********************************************************************************************************************
string SplitAbundCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta")        {  pattern = "[filename],[tag],[tag2],fasta";            }
        else if (type == "list")    {   pattern = "[filename],[tag],[tag2],list";            }
        else if (type == "name")    {   pattern = "[filename],[tag],names-[filename],[tag],[tag2],names";           }
        else if (type == "count")   {   pattern = "[filename],[tag],[tag2],count_table-[filename],[tag],count_table";     }
        else if (type == "group")   {   pattern = "[filename],[tag],[tag2],groups";          }
        else if (type == "accnos")  {   pattern = "[filename],[tag],[tag2],accnos";          }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "SplitAbundCommand", "getOutputPattern");
        exit(1);
    }
}

//**********************************************************************************************************************
SplitAbundCommand::SplitAbundCommand(string option) : Command()  {
	try {
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		else {
			OptionParser parser(option, setParameters());
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			

			//check for required parameters
			listfile = validParameter.validFile(parameters, "list");
			if (listfile == "not open") { abort = true; }
			else if (listfile == "not found") { listfile = ""; }
			else{ inputFile = listfile; current->setListFile(listfile); }	
			
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { abort = true; }
			else if (namefile == "not found") { namefile = ""; }	
			else{ inputFile = namefile; current->setNameFile(namefile); }	
		
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not open") { abort = true; }
			else if (fastafile == "not found") { fastafile = ""; }
            else { current->setFastaFile(fastafile); }
			
			groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not open") {  groupfile = ""; abort = true; }	
			else if (groupfile == "not found") { groupfile = ""; }
			else { current->setGroupFile(groupfile); }
			
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { countfile = ""; abort = true; }
			else if (countfile == "not found") { countfile = "";  }	
			else {
                current->setCountFile(countfile); 
                ct.readTable(countfile, true, false);
            }
            
            if ((namefile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: name or count.\n");  abort = true;
            }
			
            if ((groupfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: group or count.\n");  abort=true;
            }
            
			//do you have all files needed
			if ((listfile == "") && (namefile == "") && (countfile == "") && (fastafile == "")) {
				namefile = current->getNameFile(); 
				if (namefile != "") { m->mothurOut("Using " + namefile + " as input file for the name parameter.\n");  }
				else { 				
					listfile = current->getListFile(); 
					if (listfile != "") { m->mothurOut("Using " + listfile + " as input file for the list parameter.\n");  }
					else { 	
                        countfile  = current->getCountFile(); 
                        if (countfile != "") { m->mothurOut("Using " + countfile + " as input file for the count parameter.\n");  }
                        else {
                            fastafile = current->getFastaFile();
                        if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter.\n");  }
                        else {     m->mothurOut("You have no current fastafile file, and a fasta, list, name or count file is required.\n"); abort = true; } }
                    }
				}
			}
            
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = "";  allLines = true; }
			else { 
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
				else { allLines = true;  }
			}
			
			string temp = validParameter.valid(parameters, "accnos");		if (temp == "not found") { temp = "F"; }
			accnos = util.isTrue(temp); 
			
			temp = validParameter.valid(parameters, "cutoff");				if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, cutoff);

			if (cutoff == 0) {  m->mothurOut("[ERROR]: You must provide a cutoff to qualify what is abundant for the split.abund command. \n");  abort = true;  }
		}

	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "SplitAbundCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
SplitAbundCommand::~SplitAbundCommand(){}
//**********************************************************************************************************************
int SplitAbundCommand::execute(){
	try {
	
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		if (listfile != "")         { splitList();  }
		else if (namefile != "")    { splitNames(); }
        else if (countfile != "")   { splitCount(); }
		
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
		
		itTypes = outputTypes.find("accnos");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setAccnosFile(currentName); }
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
		m->errorOut(e, "SplitAbundCommand", "execute");
		exit(1);
	}
}
/**********************************************************************************************************************/
int SplitAbundCommand::splitList() {
    try {
        InputData input(listfile, "list", nullVector);
        set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        
        ListVector* list = util.getNextList(input, allLines, userLabels, processedLabels, lastLabel);
        
        if (cutoff < 1) { //percentage instead of raw count
            int total = list->getNumSeqs();
            if (countfile != "") { total = ct.getNumSeqs(); }
            else if (namefile != "") { total = util.scanNames(namefile); }
            float percentage = cutoff;
            cutoff = int(percentage * total);
            
            m->mothurOut("\nSetting cutoff to " + toString(cutoff) + "\n");
        }
        
        if (m->getControl_pressed()) { delete list; return 0; }
        
        while (list != nullptr) {
            
            if (m->getControl_pressed()) { delete list; break; }
            
            process(list); delete list;
    
            list = util.getNextList(input, allLines, userLabels, processedLabels, lastLabel);
        }
        
        if (m->getControl_pressed()) { for (int i = 0; i < outputNames.size(); i++) {    util.mothurRemove(outputNames[i]); }    return 0;    }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "SplitAbundCommand", "splitList");
        exit(1);
    }
}
/**********************************************************************************************************************/
int SplitAbundCommand::process(ListVector* thisList) {
	try {
        set<string> rareNames;
		set<string> abundNames;
        set<string> abundOTUs, rareOTUs;
		
		//get rareNames and abundNames
		for (int i = 0; i < thisList->getNumBins(); i++) {
			if (m->getControl_pressed()) { return 0; }
			
			string bin = thisList->get(i);
						
			vector<string> names;
			util.splitAtComma(bin, names);  //parses bin into individual sequence names
			int size = names.size();
            
            //if countfile is not blank we assume the list file is unique, otherwise we assume it includes all seqs
            if (countfile != "") {
                size = 0;
                for (int j = 0; j < names.size(); j++) {  size += ct.getNumSeqs(names[j]); }
            }
            
			if (size <= cutoff) {
				for (int j = 0; j < names.size(); j++) {  rareNames.insert(names[j]);  }
                rareOTUs.insert(thisList->getOTUName(i));
			}else{
				for (int j = 0; j < names.size(); j++) {  abundNames.insert(names[j]);  }
                abundOTUs.insert(thisList->getOTUName(i));
			}
		}//end for
        
        string tag = thisList->getLabel();
        vector<string> accnosOTUs = writeAccnos(tag+"_OTUS", rareOTUs, abundOTUs); //return rare, abund accnos files
        
        map<string, string> variables;
        string thisOutputDir = outputdir;
        if (outputdir == "") { thisOutputDir = util.hasPath(listfile); }
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(listfile));
        variables["[tag]"] = tag;
        variables["[tag2]"] = "rare";
        string rareList = getOutputFileName("list",variables);
        variables["[tag2]"] = "abund";
        string abundList = getOutputFileName("list",variables);
        
        string inputString = "accnos=" + accnosOTUs[0] + ", list=" + listfile + ", label=" + tag; //get rare
         m->mothurOut("/******************************************/\n");
         m->mothurOut("Running command: get.otus(" + inputString + ")\n");
         
         Command* getOTUSCommand = new GetOtusCommand(inputString);
         getOTUSCommand->execute();
         
         map<string, vector<string> > filenames = getOTUSCommand->getOutputFiles();
         
         delete getOTUSCommand;
         
         util.renameFile(filenames["list"][0], rareList);
        outputNames.push_back(rareList); outputTypes["list"].push_back(rareList);
         
         m->mothurOut("/******************************************/\nDone.\n");
         inputString = "accnos=" + accnosOTUs[1] + ", list=" + listfile+ ", label=" + tag;
         
         m->mothurOut("/******************************************/\n");
         m->mothurOut("Running command: get.otus(" + inputString + ")\n");
         
         getOTUSCommand = new GetOtusCommand(inputString);
         getOTUSCommand->execute();
         
         filenames = getOTUSCommand->getOutputFiles();
         
         delete getOTUSCommand;
         
         util.renameFile(filenames["list"][0], abundList);
         outputNames.push_back(abundList); outputTypes["list"].push_back(abundList);
         m->mothurOut("/******************************************/\nDone.\n");
         
        
        string rareCount, abundCount, rareName, abundName, rareFasta, abundFasta, rareGroup, abundGroup;
        string inputString2 = "";
        if (countfile != "") {
            if (outputdir == "") { thisOutputDir = util.hasPath(countfile); }
            variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(countfile));
            variables["[tag]"] = tag;
            variables["[tag2]"] = "rare";
            rareCount = getOutputFileName("count",variables);
            variables["[tag2]"] = "abund";
            abundCount = getOutputFileName("count",variables);
            inputString2 += ", count=" + countfile;
        }else if (groupfile != "") {
            if (outputdir == "") { thisOutputDir = util.hasPath(groupfile); }
            variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(groupfile));
            variables["[tag]"] = tag;
            variables["[tag2]"] = "rare";
            rareGroup = getOutputFileName("group",variables);
            variables["[tag2]"] = "abund";
            abundGroup = getOutputFileName("group",variables);
            inputString2 += ", group=" + groupfile;
        }
        
        if (fastafile != "")   {
            if (outputdir == "") { thisOutputDir = util.hasPath(fastafile); }
            variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(fastafile));
            variables["[tag]"] = tag;
            variables["[tag2]"] = "rare";
            rareFasta = getOutputFileName("fasta",variables);
            variables["[tag2]"] = "abund";
            abundFasta = getOutputFileName("fasta",variables);
            inputString2 += ", fasta=" + fastafile;
        }
        
        if (namefile != "")   {
            if (outputdir == "") { thisOutputDir = util.hasPath(namefile); }
            variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(namefile));
            variables["[tag]"] = tag;
            variables["[tag2]"] = "rare";
            rareName = getOutputFileName("name",variables);
            variables["[tag2]"] = "abund";
            abundName = getOutputFileName("name",variables);
            inputString2 += ", name=" + namefile;
        }
        
        if (inputString2 != "") {
            vector<string> accnosNames = writeAccnos(tag, rareNames, abundNames); //return rare, abund accnos files
            
            inputString = "dups=t, accnos=" + accnosNames[0]; //get rare
            
            m->mothurOut("/******************************************/\n");
            m->mothurOut("Running command: get.seqs(" + inputString + inputString2 + ")\n");
            
            Command* getCommand = new GetSeqsCommand(inputString+ inputString2);
            getCommand->execute();
            
            map<string, vector<string> > filenames = getCommand->getOutputFiles();
            
            delete getCommand;
            
            if (countfile != "")        {
                util.renameFile(filenames["count"][0], rareCount);
                outputNames.push_back(rareCount); outputTypes["count"].push_back(rareCount);
            }
            else if (groupfile != "")   {
                util.renameFile(filenames["group"][0], rareGroup);
                outputNames.push_back(rareGroup); outputTypes["group"].push_back(rareGroup);
            }
            if (fastafile != "")        {
                util.renameFile(filenames["fasta"][0], rareFasta);
                outputNames.push_back(rareFasta); outputTypes["fasta"].push_back(rareFasta);
            }
            if (namefile != "")         {
                util.renameFile(filenames["name"][0], rareName);
                outputNames.push_back(rareName); outputTypes["name"].push_back(rareName);
            }
            
            m->mothurOut("/******************************************/\nDone.\n");
            inputString = "dups=t, accnos=" + accnosNames[1]; //get rare
            
            m->mothurOut("/******************************************/\n");
            m->mothurOut("Running command: get.seqs(" + inputString + inputString2 + ")\n");
            
            getCommand = new GetSeqsCommand(inputString+ inputString2);
            getCommand->execute();
            
            filenames = getCommand->getOutputFiles();
            
            delete getCommand;
            
            if (countfile != "")        { util.renameFile(filenames["count"][0], abundCount);  outputNames.push_back(abundCount); outputTypes["count"].push_back(abundCount);  }
            else if (groupfile != "")   { util.renameFile(filenames["group"][0], abundGroup);
                outputNames.push_back(abundGroup); outputTypes["group"].push_back(abundGroup);
            }
            if (fastafile != "")        { util.renameFile(filenames["fasta"][0], abundFasta); outputNames.push_back(abundFasta); outputTypes["fasta"].push_back(abundFasta);   }
            if (namefile != "")         { util.renameFile(filenames["name"][0], abundName);   outputNames.push_back(abundName); outputTypes["name"].push_back(abundName);   }
            
            m->mothurOut("/******************************************/\nDone.\n");
        }
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "process");
		exit(1);
	}
}
/**********************************************************************************************************************/
int SplitAbundCommand::splitCount() { //countfile
	try {
        inputFile = countfile;
		set<string> rareNames;
		set<string> abundNames;
        
        if (cutoff < 1) { //cutoff is a percentage rather than a explicit size
            float percentage = cutoff;
            int totalSeqs = ct.getNumSeqs();
            
            cutoff = int(totalSeqs * percentage);
            
            m->mothurOut("\nSetting cutoff to " + toString(cutoff) + "\n");
        }
        
		vector<string> allNames = ct.getNamesOfSeqs();
        for (int i = 0; i < allNames.size(); i++) {
            
            if (m->getControl_pressed()) { return 0; }
            
            int size = ct.getNumSeqs(allNames[i]);
            
			if (size <= cutoff) {
				rareNames.insert(allNames[i]); 
			}else{
				abundNames.insert(allNames[i]); 
			}
		}
        
        vector<string> accnosNames = writeAccnos("", rareNames, abundNames); //return rare, abund accnos files
        
        map<string, string> variables;
        string thisOutputDir = outputdir;
        if (outputdir == "") { thisOutputDir = util.hasPath(fastafile); }
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(fastafile));
        variables["[tag]"] = "";
        variables["[tag2]"] = "rare";
        string rareFasta = getOutputFileName("fasta",variables);
        variables["[tag2]"] = "abund";
        string abundFasta = getOutputFileName("fasta",variables);
        
        if (outputdir == "") { thisOutputDir = util.hasPath(countfile); }
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(countfile));
        variables["[tag]"] = "";
        variables["[tag2]"] = "rare";
        string rareCount = getOutputFileName("count",variables);
        variables["[tag2]"] = "abund";
        string abundCount = getOutputFileName("count",variables);
        
        
        string inputString = "dups=t, accnos=" + accnosNames[0] + ", count=" + countfile; //get rare
        if (fastafile != "")   { inputString += ", fasta=" + fastafile;      }
        
        m->mothurOut("/******************************************/\n");
        m->mothurOut("Running command: get.seqs(" + inputString + ")\n");
        
        Command* getCommand = new GetSeqsCommand(inputString);
        getCommand->execute();
        
        map<string, vector<string> > filenames = getCommand->getOutputFiles();
        
        delete getCommand;
        
        util.renameFile(filenames["count"][0], rareCount);
        outputNames.push_back(rareCount); outputTypes["count"].push_back(rareCount);
        if (fastafile != "")   { util.renameFile(filenames["fasta"][0], rareFasta); outputNames.push_back(rareFasta); outputTypes["fasta"].push_back(rareFasta); }
        
        m->mothurOut("/******************************************/\nDone.\n");
        inputString = "dups=t, accnos=" + accnosNames[1] + ", count=" + countfile; //get rare
        if (fastafile != "")   { inputString += ", fasta=" + fastafile;      }
        
        m->mothurOut("/******************************************/\n");
        m->mothurOut("Running command: get.seqs(" + inputString + ")\n");
        
        getCommand = new GetSeqsCommand(inputString);
        getCommand->execute();
        
        filenames = getCommand->getOutputFiles();
        
        delete getCommand;
        
        util.renameFile(filenames["count"][0], abundCount);
        outputNames.push_back(abundCount); outputTypes["count"].push_back(abundCount);
        if (fastafile != "")   { util.renameFile(filenames["fasta"][0], abundFasta); outputNames.push_back(abundFasta); outputTypes["fasta"].push_back(abundFasta);  }
        
        m->mothurOut("/******************************************/\nDone.\n");
		
		return 0;  
	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "splitCount");
		exit(1);
	}
}
/**********************************************************************************************************************/
int SplitAbundCommand::splitNames() { //namefile
	try {
		set<string> rareNames;
		set<string> abundNames;
        
        if (cutoff < 1) { //percentage used, find cutoff value
            float percentage = cutoff;
            int totalSeqs = util.scanNames(namefile);
            
            cutoff = int(totalSeqs * percentage);
            
            m->mothurOut("\nSetting cutoff to " + toString(cutoff) + "\n");
        }
			
		//open input file
		ifstream in; util.openInputFile(namefile, in);
		
		while (!in.eof()) {
			if (m->getControl_pressed()) { break; }
			
			string firstCol, secondCol;
            in >> firstCol; util.gobble(in); in >> secondCol; util.gobble(in);
			
			int size = util.getNumNames(secondCol);
				
			if (size <= cutoff) {
				rareNames.insert(firstCol); 
			}else{
				abundNames.insert(firstCol); 
			}
		}
		in.close();
        
        vector<string> accnosNames = writeAccnos("", rareNames, abundNames); //return rare, abund accnos files
    
        map<string, string> variables;
        string thisOutputDir = outputdir;
        if (outputdir == "") { thisOutputDir = util.hasPath(groupfile); }
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(groupfile));
        variables["[tag]"] = "";
        variables["[tag2]"] = "rare";
        string rareGroup = getOutputFileName("group",variables);
        variables["[tag2]"] = "abund";
        string abundGroup = getOutputFileName("group",variables);
        
        if (outputdir == "") { thisOutputDir = util.hasPath(fastafile); }
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(fastafile));
        variables["[tag]"] = "";
        variables["[tag2]"] = "rare";
        string rareFasta = getOutputFileName("fasta",variables);
        variables["[tag2]"] = "abund";
        string abundFasta = getOutputFileName("fasta",variables);
        
        if (outputdir == "") { thisOutputDir = util.hasPath(namefile); }
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(namefile));
        variables["[tag]"] = "";
        variables["[tag2]"] = "rare";
        string rareName = getOutputFileName("name",variables);
        variables["[tag2]"] = "abund";
        string abundName = getOutputFileName("name",variables);
        
        
        string inputString = "dups=t, accnos=" + accnosNames[0] + ", name=" + namefile; //get rare
        if (groupfile != "")   { inputString += ", group=" + groupfile;      }
        if (fastafile != "")   { inputString += ", fasta=" + fastafile;      }
        
        m->mothurOut("/******************************************/\n");
        m->mothurOut("Running command: get.seqs(" + inputString + ")\n");
        
        Command* getCommand = new GetSeqsCommand(inputString);
        getCommand->execute();
        
        map<string, vector<string> > filenames = getCommand->getOutputFiles();
        
        delete getCommand;
        
        util.renameFile(filenames["name"][0], rareName);
        outputNames.push_back(rareName); outputTypes["name"].push_back(rareName);
        if (groupfile != "")   { util.renameFile(filenames["group"][0], rareGroup); outputNames.push_back(rareGroup); outputTypes["group"].push_back(rareGroup); }
        if (fastafile != "")   { util.renameFile(filenames["fasta"][0], rareFasta); outputNames.push_back(rareFasta); outputTypes["fasta"].push_back(rareFasta); }
        
        m->mothurOut("/******************************************/\nDone.\n");
        
        inputString = "dups=t, accnos=" + accnosNames[1] + ", name=" + namefile; //get rare
        if (groupfile != "")   { inputString += ", group=" + groupfile;      }
        if (fastafile != "")   { inputString += ", fasta=" + fastafile;      }
        
        m->mothurOut("/******************************************/\n");
        m->mothurOut("Running command: get.seqs(" + inputString + ")\n");
        
        getCommand = new GetSeqsCommand(inputString);
        getCommand->execute();
        
        filenames = getCommand->getOutputFiles();
        
        delete getCommand;
        
        util.renameFile(filenames["name"][0], abundName);
        outputNames.push_back(abundName); outputTypes["name"].push_back(abundName);
        if (groupfile != "")   { util.renameFile(filenames["group"][0], abundGroup); outputNames.push_back(abundGroup); outputTypes["group"].push_back(abundGroup); }
        if (fastafile != "")   { util.renameFile(filenames["fasta"][0], abundFasta); outputNames.push_back(abundFasta); outputTypes["fasta"].push_back(abundFasta); }
        
        m->mothurOut("/******************************************/\nDone.\n");
				
		return 0;

	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "splitNames");
		exit(1);
	}
}
/**********************************************************************************************************************/
//just write the unique names - if a namesfile is given
vector<string> SplitAbundCommand::writeAccnos(string tag, set<string> rareNames, set<string> abundNames) {
	try {
        vector<string> outputAccnosFiles;
        map<string, string> variables;
        string thisOutputDir = outputdir;
        if (outputdir == "") { thisOutputDir = util.hasPath(inputFile); }
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(inputFile));
        variables["[tag]"] = tag;
        variables["[tag2]"] = "rare";
        string rare = getOutputFileName("accnos",variables);
        outputAccnosFiles.push_back(rare);
        
        ofstream aout, rout;
        util.openOutputFile(rare, rout);
        outputNames.push_back(rare); outputTypes["accnos"].push_back(rare);
        
        for (set<string>::iterator itRare = rareNames.begin(); itRare != rareNames.end(); itRare++) { rout << (*itRare) << endl; }
        rout.close();
        
        variables["[tag2]"] = "abund";
        string abund = getOutputFileName("accnos",variables);
        util.openOutputFile(abund, aout);
        outputNames.push_back(abund); outputTypes["accnos"].push_back(abund);
        outputAccnosFiles.push_back(abund);
        
        for (set<string>::iterator itAbund = abundNames.begin(); itAbund != abundNames.end(); itAbund++) { aout << (*itAbund) << endl; }
        aout.close();

        return outputAccnosFiles;

	}
	catch(exception& e) {
		m->errorOut(e, "SplitAbundCommand", "writeAccnos");
		exit(1);
	}
}
/**********************************************************************************************************************/
