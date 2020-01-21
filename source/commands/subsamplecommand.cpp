/*
 *  subsamplecommand.cpp
 *  Mothur
 *
 *  Created by westcott on 10/27/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "subsamplecommand.h"

#include "deconvolutecommand.h"
#include "getseqscommand.h"
#include "subsample.h"

//**********************************************************************************************************************
vector<string> SubSampleCommand::setParameters(){	
	try {		
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "FLSSR", "none","fasta",false,false,true); parameters.push_back(pfasta);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","name",false,false,true); parameters.push_back(pname);
        CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "none", "none","taxonomy",false,false,true); parameters.push_back(ptaxonomy);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","count",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","group",false,false,true); parameters.push_back(pgroup);
		CommandParameter plist("list", "InputTypes", "", "", "none", "FLSSR", "none","list",false,false,true); parameters.push_back(plist);
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "FLSSR", "none","shared",false,false,true); parameters.push_back(pshared);
		CommandParameter prabund("rabund", "InputTypes", "", "", "none", "FLSSR", "none","rabund",false,false); parameters.push_back(prabund);
		CommandParameter psabund("sabund", "InputTypes", "", "", "none", "FLSSR", "none","sabund",false,false); parameters.push_back(psabund);
        CommandParameter ptree("tree", "InputTypes", "", "", "none", "FLSSR", "none","tree",false,false); parameters.push_back(ptree);
        CommandParameter pconstaxonomy("constaxonomy", "InputTypes", "", "", "none", "none", "none","constaxonomy",false,false, true); parameters.push_back(pconstaxonomy);
		CommandParameter plabel("label", "String", "", "", "", "", "","",false,false); parameters.push_back(plabel);
		CommandParameter pgroups("groups", "String", "", "", "", "", "","",false,false); parameters.push_back(pgroups);
		CommandParameter psize("size", "Number", "", "0", "", "", "","",false,false,true); parameters.push_back(psize);
		CommandParameter ppersample("persample", "Boolean", "", "F", "", "", "","",false,false,true); parameters.push_back(ppersample);
        CommandParameter pwithreplacement("withreplacement", "Boolean", "", "F", "", "", "","",false,false,true); parameters.push_back(pwithreplacement);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SubSampleCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The sub.sample command is designed to be used as a way to normalize your data, or create a smaller set from your original set.\n";
		helpString += "The sub.sample command parameters are " + getCommandParameters() + ".  You must provide a fasta, list, sabund, rabund or shared file as an input file.\n";
		helpString += "The namefile is only used with the fasta file, not with the listfile, because the list file should contain all sequences.\n";
		helpString += "The groups parameter allows you to specify which of the groups in your groupfile you would like included. The group names are separated by dashes.\n";
		helpString += "The label parameter allows you to select what distance levels you would like, and are also separated by dashes.\n";
		helpString += "The size parameter allows you indicate the size of your subsample.\n";
		helpString += "The persample parameter allows you indicate you want to select subsample of the same size from each of your groups, default=false. It is only used with the list and fasta files if a groupfile is given.\n";
		helpString += "persample=false will select a random set of sequences of the size you select, but the number of seqs from each group may differ.\n";
		helpString += "The size parameter is not set: with shared file size=number of seqs in smallest sample, with all other files if a groupfile is given and persample=true, then size=number of seqs in smallest sample, otherwise size=10% of number of seqs.\n";
        helpString += "The withreplacement parameter allows you to indicate you want to subsample your data allowing for the same read to be included multiple times. Default=f. \n";
		helpString += "The sub.sample command should be in the following format: sub.sample(list=yourListFile, group=yourGroupFile, groups=yourGroups, label=yourLabels).\n";
		helpString += "Example sub.sample(list=abrecovery.fn.list, group=abrecovery.groups, groups=B-C, size=20).\n";
		helpString += "The default value for groups is all the groups in your groupfile, and all labels in your inputfile will be used.\n";
		helpString += "The sub.sample command outputs a .subsample file.\n";
		
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string SubSampleCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta")            {   pattern = "[filename],subsample,[extension]";    }
        else if (type == "sabund")      {   pattern = "[filename],subsample,[extension]";    }
        else if (type == "name")        {   pattern = "[filename],subsample,[extension]";    }
        else if (type == "group")       {   pattern = "[filename],subsample,[extension]";    }
        else if (type == "count")       {   pattern = "[filename],subsample,[extension]";    }
        else if (type == "tree")        {   pattern = "[filename],subsample,[extension]";    }
        else if (type == "list")        {   pattern = "[filename],[distance],subsample,[extension]";    }
        else if (type == "taxonomy")    {   pattern = "[filename],subsample,[extension]";    }
        else if (type == "constaxonomy"){   pattern = "[filename],subsample,[extension]";    }
        else if (type == "shared")      {   pattern = "[filename],[distance],subsample,[extension]";    }
        else if (type == "rabund")      {   pattern = "[filename],subsample,[extension]";    }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "SubSampleCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
SubSampleCommand::SubSampleCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["shared"] = tempOutNames;
		outputTypes["list"] = tempOutNames;
		outputTypes["rabund"] = tempOutNames;
		outputTypes["sabund"] = tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["name"] = tempOutNames;
		outputTypes["group"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
        outputTypes["taxonomy"] = tempOutNames;
        outputTypes["constaxonomy"] = tempOutNames;
        outputTypes["tree"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "GetRelAbundCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
SubSampleCommand::SubSampleCommand(string option) {
	try {
		abort = false; calledHelp = false;   
		allLines = true;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			
			//check to make sure all parameters are valid for command
			map<string,string>::iterator it;
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["shared"] = tempOutNames;
			outputTypes["list"] = tempOutNames;
			outputTypes["rabund"] = tempOutNames;
			outputTypes["sabund"] = tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["name"] = tempOutNames;
			outputTypes["group"] = tempOutNames;
            outputTypes["count"] = tempOutNames;
            outputTypes["taxonomy"] = tempOutNames;
            outputTypes["constaxonomy"] = tempOutNames;
            outputTypes["tree"] = tempOutNames;
					
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";	}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
				
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
				
				it = parameters.find("sabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["sabund"] = inputDir + it->second;		}
				}
				
				it = parameters.find("rabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["rabund"] = inputDir + it->second;		}
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
                
                it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
				}
                
                it = parameters.find("constaxonomy");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["constaxonomy"] = inputDir + it->second;		}
                }

                
                it = parameters.find("tree");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["tree"] = inputDir + it->second;		}
                }

			}
			
			//check for required parameters
			listfile = validParameter.validFile(parameters, "list");
			if (listfile == "not open") { listfile = ""; abort = true; }
			else if (listfile == "not found") { listfile = ""; }
			else { current->setListFile(listfile); }
			
			sabundfile = validParameter.validFile(parameters, "sabund");
			if (sabundfile == "not open") { sabundfile = ""; abort = true; }	
			else if (sabundfile == "not found") { sabundfile = ""; }
			else { current->setSabundFile(sabundfile); }
			
			rabundfile = validParameter.validFile(parameters, "rabund");
			if (rabundfile == "not open") { rabundfile = ""; abort = true; }	
			else if (rabundfile == "not found") { rabundfile = ""; }
			else { current->setRabundFile(rabundfile); }
			
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not open") { fastafile = ""; abort = true; }	
			else if (fastafile == "not found") { fastafile = ""; }
			else { current->setFastaFile(fastafile); }
			
			sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { sharedfile = ""; abort = true; }	
			else if (sharedfile == "not found") { sharedfile = ""; }
			else { current->setSharedFile(sharedfile); }
			
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { namefile = ""; abort = true; }	
			else if (namefile == "not found") { namefile = ""; }
			else { current->setNameFile(namefile); }
			
			groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not open") { groupfile = ""; abort = true; }	
			else if (groupfile == "not found") { groupfile = ""; }
			else { current->setGroupFile(groupfile); }
            
            taxonomyfile = validParameter.validFile(parameters, "taxonomy");
			if (taxonomyfile == "not open") { taxonomyfile = ""; abort = true; }
			else if (taxonomyfile == "not found") { taxonomyfile = ""; }
			else { current->setTaxonomyFile(taxonomyfile); }
            
            constaxonomyfile = validParameter.validFile(parameters, "constaxonomy");
            if (constaxonomyfile == "not open") { constaxonomyfile = ""; abort = true; }
            else if (constaxonomyfile == "not found") { constaxonomyfile = ""; }
            else { current->setConsTaxonomyFile(constaxonomyfile); }
            
            treefile = validParameter.validFile(parameters, "tree");
            if (treefile == "not open") { treefile = ""; abort = true; }
            else if (treefile == "not found") { treefile = ""; }
            else { current->setTreeFile(treefile); }
			
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
            
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			label = validParameter.valid(parameters, "label");			
			if (label == "not found") { label = ""; }
			else { 
				if(label != "all") {  util.splitAtDash(label, labels);  allLines = false;  }
				else { allLines = true;  }
			}
			
			groups = validParameter.valid(parameters, "groups");			
			if (groups == "not found") { groups = ""; pickedGroups = false; }
			else { 
				pickedGroups = true;
				util.splitAtDash(groups, Groups);
                if (Groups.size() != 0) { if (Groups[0]== "all") { Groups.clear(); } }
			}
			
			string temp = validParameter.valid(parameters, "size");		if (temp == "not found"){	temp = "0";		}
			util.mothurConvert(temp, size);  
			
			temp = validParameter.valid(parameters, "persample");		if (temp == "not found"){	temp = "f";		}
			persample = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "withreplacement");		if (temp == "not found"){	temp = "f";		}
            withReplacement = util.isTrue(temp);
			
			if ((groupfile == "") && (countfile == "")) { persample = false; }
            if (countfile != "") {
                if (!ct.hasGroupInfo()) { 
                    persample = false; 
                    if (pickedGroups) { m->mothurOut("You cannot pick groups without group info in your count file.\n");  abort = true; }
                }
            }
			
			if ((namefile != "") && ((fastafile == "") && (taxonomyfile == "") && (treefile == ""))) { m->mothurOut("You may only use a name file with a fasta file, tree file or taxonomy file.\n");  abort = true; }
            
            if ((taxonomyfile != "") && ((fastafile == "") && (listfile == ""))) { m->mothurOut("You may only use a taxonomy file with a fasta file or list file.\n");  abort = true; }
            if ((constaxonomyfile != "") && ((sharedfile == "") && (listfile == ""))) { m->mothurOut("You may only use a constaxonomy file with a shared file or list file.\n");  abort = true; }
            
			
			if ((fastafile == "") && (listfile == "") && (sabundfile == "") && (rabundfile == "") && (sharedfile == "") && (treefile == "")) {
				m->mothurOut("You must provide a fasta, list, sabund, rabund, shared or tree file as an input file.\n"); abort = true; }
			
			if (pickedGroups && ((groupfile == "") && (sharedfile == "") && (countfile == ""))) { 
				m->mothurOut("You cannot pick groups without a valid group, count or shared file.\n");  abort = true; }
			
			if (((groupfile != "") || (countfile != "")) && ((fastafile == "") && (listfile == "") && (treefile == ""))) {
				m->mothurOut("Group or count files are only valid with list file, fasta file or tree file.\n"); abort = true; }
			
			if (((groupfile != "") || (countfile != "")) && ((fastafile != "") && (listfile != ""))) { 
				m->mothurOut("A new group or count file can only be made from the subsample of a list file or fasta file, not both. Please correct.\n"); abort = true; }
			
            if (countfile == "") {
                if ((fastafile != "") && (namefile == "")) {
                    vector<string> files; files.push_back(fastafile);
                    if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
                }
            }
		}

	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "SubSampleCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int SubSampleCommand::execute(){
	try {
	
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		if (sharedfile != "")       {   getSubSampleShared();	}
		else if (listfile != "")	{   getSubSampleList();		} //can only use with replacement if a count file is provided
		else if (rabundfile != "")	{   getSubSampleRabund();	}
		else if (sabundfile != "")	{   getSubSampleSabund();	}
		else if (fastafile != "")	{   getSubSampleFasta();	} //can only use with replacement if a count file is provided
		else if (treefile != "")	{   getSubSampleTree();     }
        
        if (m->getControl_pressed()) {  for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);  } return 0; }
			
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
		
		itTypes = outputTypes.find("shared");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setSharedFile(currentName); }
		}
		
		itTypes = outputTypes.find("rabund");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setRabundFile(currentName); }
		}
		
		itTypes = outputTypes.find("sabund");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setSabundFile(currentName); }
		}
        
        itTypes = outputTypes.find("count");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setCountFile(currentName); }
		}
		
        itTypes = outputTypes.find("taxonomy");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setTaxonomyFile(currentName); }
		}
        
        itTypes = outputTypes.find("constaxonomy");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setConsTaxonomyFile(currentName); }
        }
        
        itTypes = outputTypes.find("tree");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setTreeFile(currentName); }
        }
		
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int SubSampleCommand::getSubSampleTree() {
    try {
        TreeReader* reader = NULL;
        string cinputfile = countfile;
        if (countfile == "") {
            reader = new TreeReader(treefile, groupfile, namefile);
            cinputfile = namefile;
        }else {
            reader = new TreeReader(treefile, countfile);
        }
        
        vector<Tree*> T; T = reader->getTrees();   //user trees
        CountTable* ct; ct = T[0]->getCountTable();
        vector<string> Treenames = T[0]->getTreeNames();
        delete reader;
        
        if (size == 0) { //user has not set size, set size = smallest samples size
            size = ct->getNumSeqsSmallestGroup();
        }else {
            ct->removeGroup(size);
        }
        
        if (!pickedGroups) {  Groups = ct->getNamesOfGroups();  }
        if (Groups.size() == 0) {  m->mothurOut("The size you selected is too large, skipping tree file.\n");   return 0; }
        
        m->mothurOut("Sampling " + toString(size) + " from each group.\n");
        
        //copy to preserve old one - would do this in subsample but memory cleanup becomes messy.
        CountTable* newCt = new CountTable();
        
        SubSample sample;
        Tree* subSampleTree; //sets unwanted seqs to doNotIncludeMe, prune below
        if (withReplacement)    { subSampleTree = sample.getSampleWithReplacement(T[0], ct, newCt, size, Groups);       }
        else                    { subSampleTree = sample.getSample(T[0], ct, newCt, size, Groups);                      }
         
        if (m->getControl_pressed()) { delete newCt; delete subSampleTree; return 0; }
        
        vector<string> newTreeNames = newCt->getNamesOfSeqs(Groups); //without "doNotIncludeMe"
        Tree* outputTree = new Tree(newTreeNames.size(), ct, Treenames);  //create ouptut tree - respecting pickedGroups
        
        outputTree->getSubTree(subSampleTree, newTreeNames); //builds tree with only seqs present in newCt
        outputTree->assembleTree();
        newCt->removeGroup("doNotIncludeMe");
        
        string treeOutputDir = outputDir;
        if (outputDir == "") {  treeOutputDir += util.hasPath(treefile);  }
        map<string, string> variables;
        variables["[filename]"] = treeOutputDir + util.getRootName(util.getSimpleName(treefile));
        variables["[extension]"] = util.getExtension(treefile);
        string treeOutputFileName = getOutputFileName("tree", variables);
        outputTypes["tree"].push_back(treeOutputFileName);  outputNames.push_back(treeOutputFileName);
        
        ofstream outTree; util.openOutputFile(treeOutputFileName, outTree);
        outputTree->print(outTree, "both"); outTree.close();
        for (int i = 0; i < T.size(); i++) { delete T[i]; } delete subSampleTree; delete outputTree;
        
        string countOutputDir = outputDir;
        if (outputDir == "") {  countOutputDir += util.hasPath(cinputfile);  }
        variables["[filename]"] = countOutputDir + util.getRootName(util.getSimpleName(cinputfile));
        variables["[extension]"] = ".count_table";
        string countOutputFile = getOutputFileName("count",variables);
        outputNames.push_back(countOutputFile); outputTypes["count"].push_back(countOutputFile);
        
        newCt->printTable(countOutputFile); delete newCt;
        delete ct;
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "SubSampleCommand", "getSubSampleTree");
        exit(1);
    }
}
//**********************************************************************************************************************
int SubSampleCommand::getSubSampleFasta() {
	try {
        vector<string> names;
		if (namefile != "") { names = readNames(); }	//fills names with all names in namefile.
		else { names = getNames(); }//no name file, so get list of names to pick from
		
		GroupMap groupMap;
		if (groupfile != "") {
			groupMap.readMap(groupfile);
			
            if (Groups.size() == 0) { Groups = groupMap.getNamesOfGroups(); }
			//file mismatch quit
			if (names.size() != groupMap.getNumSeqs()) { 
				m->mothurOut("[ERROR]: your fasta file contains " + toString(names.size()) + " sequences, and your groupfile contains " + toString(groupMap.getNumSeqs()) + ", please correct.\n"); return 0;
			}			
		}else if (countfile != "") {
            if (ct.hasGroupInfo()) { if (Groups.size() == 0) { Groups = ct.getNamesOfGroups(); } }
            
            //file mismatch quit
			if (names.size() != ct.getNumUniqueSeqs()) { 
				m->mothurOut("[ERROR]: your fasta file contains " + toString(names.size()) + " sequences, and your count file contains " + toString(ct.getNumUniqueSeqs()) + " unique sequences, please correct.\n"); return 0;
			}	
        }
		
		if (m->getControl_pressed()) { return 0; }
		
		//make sure that if your picked groups size is not too big
		int thisSize = 0;
        if (countfile == "") {
            thisSize = names.size();
            if (withReplacement) {  m->mothurOut("[WARNING]: To use the withreplacement option with a fasta file, you need to provide a count file, running without replacement.\n"); withReplacement = false; }
        }
        else {  thisSize = ct.getNumSeqs();  }  //all seqs not just unique
        
		if (persample) { 
			if (size == 0) { //user has not set size, set size = smallest samples size
				if (countfile == "") { size = groupMap.getNumSeqsSmallestGroup(); }
                else {  size = ct.getNumSeqsSmallestGroup();  }
            }else { //make sure size is not too large
				vector<string> newGroups;
				for (int i = 0; i < Groups.size(); i++) {
					int thisSize = 0;
                    if (countfile == "") { thisSize = groupMap.getNumSeqs(Groups[i]); }
                    else {  thisSize = ct.getGroupCount(Groups[i]);  }
					
					if (thisSize >= size) {	newGroups.push_back(Groups[i]);	}
					else {  m->mothurOut("You have selected a size that is larger than " + Groups[i] + " number of sequences, removing " + Groups[i] + ".\n");  }
				}
				Groups = newGroups;
                if (newGroups.size() == 0) {  m->mothurOut("[ERROR]: all groups removed.\n");  m->setControl_pressed(true); }
			}
			
			m->mothurOut("Sampling " + toString(size) + " from each group.\n");
		}else {
			if (pickedGroups) {
				int total = 0;
				for(int i = 0; i < Groups.size(); i++) {
                    if (countfile == "") { total += groupMap.getNumSeqs(Groups[i]); }
                    else {  total += ct.getGroupCount(Groups[i]);  }
				}
				
				if (size == 0) { size = int (total * 0.10); } //user has not set size, set size = 10% samples size
				
				if (total < size) { 
					if (size != 0) {  m->mothurOut("Your size is too large for the number of groups you selected. Adjusting to " + toString(int (total * 0.10)) + ".\n"); }
					size = int (total * 0.10);
				}
				m->mothurOut("Sampling " + toString(size) + " from " + toString(total) + ".\n");
			}
			
			if (size == 0) { //user has not set size, set size = 10% samples size
				if (countfile == "") {  size = int (names.size() * 0.10); }
                else {  size = int (ct.getNumSeqs() * 0.10); }
			}
			
            if (size > thisSize) { m->mothurOut("Your fasta file only contains " + toString(thisSize) + " sequences. Setting size to " + toString(thisSize) + ".\n");
                    size = thisSize; }
            
            if (!pickedGroups) { m->mothurOut("Sampling " + toString(size) + " from " + toString(thisSize) + ".\n");  }
		}
		util.mothurRandomShuffle(names);
		
		set<string> subset; //subset may contain names from column 2 of namefile if namefile is used. Will need to be sure to match these with unique name in fasta file
        
        if (countfile == "") { //fill subset with the names we want to sample
            SubSample sample;
            if (groupfile != "") {
                GroupMap sampledGm = sample.getSample(groupMap, size, Groups, persample);
                vector<string> sampledSeqs = sampledGm.getNamesSeqs();
                for (int i = 0; i < sampledSeqs.size(); i++) { subset.insert(sampledSeqs[i]); }
            
                string groupOutputDir = outputDir;
                if (outputDir == "") {  groupOutputDir += util.hasPath(groupfile);  }
                map<string, string> variables;
                variables["[filename]"] = groupOutputDir + util.getRootName(util.getSimpleName(groupfile));
                variables["[extension]"] = util.getExtension(groupfile);
                string groupOutputFileName = getOutputFileName("group", variables);
                outputTypes["group"].push_back(groupOutputFileName);  outputNames.push_back(groupOutputFileName);
                sampledGm.print(groupOutputFileName);
            }else {
                for (int i = 0; i < size; i++) { subset.insert(names[i]); }
            }
        }else {
            SubSample sample;
            CountTable sampledCt;
            if (withReplacement)    {  sampledCt = sample.getSampleWithReplacement(ct, size, Groups, persample);    }
            else                    {  sampledCt = sample.getSample(ct, size, Groups, persample);                   }
            
            vector<string> sampledSeqs = sampledCt.getNamesOfSeqs();
            for (int i = 0; i < sampledSeqs.size(); i++) { subset.insert(sampledSeqs[i]); }
            
            string countOutputDir = outputDir;
            if (outputDir == "") {  countOutputDir += util.hasPath(countfile);  }
            map<string, string> variables;
            variables["[filename]"] = countOutputDir + util.getRootName(util.getSimpleName(countfile));
            variables["[extension]"] = util.getExtension(countfile);
            string countOutputFileName = getOutputFileName("count", variables);
            outputTypes["count"].push_back(countOutputFileName);  outputNames.push_back(countOutputFileName);
            sampledCt.printTable(countOutputFileName);
        }
        
		
		if (subset.size() == 0) {  m->mothurOut("The size you selected is too large, skipping fasta file.\n");   return 0; }
		
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(fastafile);  }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(fastafile));
        variables["[extension]"] = util.getExtension(fastafile);
		string outputFileName = getOutputFileName("fasta", variables);
        outputTypes["fasta"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        
		ofstream out;
		util.openOutputFile(outputFileName, out);
		
		//read through fasta file outputting only the names on the subsample list
		ifstream in;
		util.openInputFile(fastafile, in);
		
		string thisname;
		int count = 0;
		map<string, vector<string> >::iterator itNameMap;
		
		while(!in.eof()){
			
			if (m->getControl_pressed()) { in.close(); out.close();  return 0; }
			
			Sequence currSeq(in);
			thisname = currSeq.getName();
			
			if (thisname != "") {
				
				//does the subset contain a sequence that this sequence represents
				itNameMap = nameMap.find(thisname);
				if (itNameMap != nameMap.end()) {
					vector<string> nameRepresents = itNameMap->second;
				
					for (int i = 0; i < nameRepresents.size(); i++){
						if (subset.count(nameRepresents[i]) != 0) {
							out << ">" << nameRepresents[i] << endl << currSeq.getAligned() << endl;
							count++;
						}
					}
				}else{ m->mothurOut("[ERROR]: " + thisname + " is not in your namefile, please correct.\n");  }
			}
			util.gobble(in);
		}
		in.close();	
		out.close();
        
		
		if (count != subset.size()) { m->mothurOut("[ERROR]: The subset selected contained " + toString(subset.size()) + " sequences, but I only found " + toString(count) + " of those in the fastafile.\n"); }
		
		if (namefile != "") {
			m->mothurOut("Deconvoluting subsampled fasta file... \n");
			map<string, string> variables; 
            variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(namefile));
            variables["[extension]"] = util.getExtension(namefile);
            string outputNameFileName = getOutputFileName("name", variables);
			//use unique.seqs to create new name and fastafile
			string inputString = "fasta=" + outputFileName;
			m->mothurOut("/******************************************/\n");
			m->mothurOut("Running command: unique.seqs(" + inputString + ")\n");
			current->setMothurCalling(true);
            
			Command* uniqueCommand = new DeconvoluteCommand(inputString);
			uniqueCommand->execute();
			
			map<string, vector<string> > filenames = uniqueCommand->getOutputFiles();
			
			delete uniqueCommand;
			current->setMothurCalling(false);
            
            util.renameFile(filenames["name"][0], outputNameFileName); 
            util.renameFile(filenames["fasta"][0], outputFileName);  
            
			outputTypes["name"].push_back(outputNameFileName);  outputNames.push_back(outputNameFileName);

			m->mothurOut("/******************************************/\n");
			
			m->mothurOut("Done.\n");
            
            if (taxonomyfile != "") {
                set<string> tempSubset;
                //get new unique names from fasta file
                //read through fasta file outputting only the names on the subsample list after deconvolute
                ifstream in2;
                util.openInputFile(outputFileName, in2);
                
                while (!in2.eof()) {
                    Sequence seq(in2); util.gobble(in2);
                    if (seq.getName() != "") {
                        tempSubset.insert(seq.getName());
                    }
                }
                in2.close();
                
                //send that list to getTax
                int tcount = getTax(tempSubset);
                if (tcount != tempSubset.size()) { m->mothurOut("[ERROR]: subsampled fasta file contains " + toString(tempSubset.size()) + " sequences, but I only found " + toString(tcount) + " in your taxonomy file, please correct.\n"); }
            }
		}else {
            if (taxonomyfile != "") {
                int tcount = getTax(subset);
                if (tcount != subset.size()) { m->mothurOut("[ERROR]: subsampled fasta file contains " + toString(subset.size()) + " sequences, but I only found " + toString(tcount) + " in your taxonomy file, please correct.\n");  }
            
            }  //should only contain uniques.
        }
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getSubSampleFasta");
		exit(1);
	}
}
//**********************************************************************************************************************
vector<string> SubSampleCommand::getNames() {
	try {
		vector<string> names;
        
		ifstream in;
		util.openInputFile(fastafile, in);
		
		string thisname;
		while(!in.eof()){
			
			if (m->getControl_pressed()) { in.close(); return names; }
			
			Sequence currSeq(in);
			thisname = currSeq.getName();
			
			if (thisname != "") {
				vector<string> temp; temp.push_back(thisname);
				nameMap[thisname] = temp;
				names.push_back(thisname);
			}
			util.gobble(in);
		}
		in.close();	
		
		return names;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getNames");
		exit(1);
	}
}	
//**********************************************************************************************************************
vector<string> SubSampleCommand::readNames() {
	try {
        vector<string> names;
        
        nameMap.clear();
        util.readNames(namefile, nameMap);
        
        //save names of all sequences
        map<string, vector<string> >::iterator it;
        for (it = nameMap.begin(); it != nameMap.end(); it++) { for (int i = 0; i < (it->second).size(); i++) { names.push_back((it->second)[i]); } }
        
		return names;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "readNames");
		exit(1);
	}
}		
//**********************************************************************************************************************
int SubSampleCommand::getSubSampleShared() {
	try {
		InputData input(sharedfile, "sharedfile", Groups);
		set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        
        SharedRAbundVectors* lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
        Groups = lookup->getNamesGroups();
        
        if (constaxonomyfile != "") { //you have a constaxonomy file and have not set the labels parameter. Either your sharedfile has one label or we only want to use one since the constaxonomy file is related to one label
            if (labels.size() == 0) {
                labels.insert(lastLabel); allLines = false;
                m->mothurOut("\n[WARNING]: The constaxonomy file represents a single label in your shared file. You did not set the label parameter, so mothur is assuming you want to use label " + lastLabel + ". If this is not correct file mismatches can occur.\n\n");
            }
        }
        
		if (size == 0) { //user has not set size, set size = smallest samples size
			size = lookup->getNumSeqsSmallestGroup();
        }else {
            lookup->removeGroups(size);
            Groups = lookup->getNamesGroups();
		}
		if (lookup->size() == 0) {  m->mothurOut("The size you selected is too large, skipping shared file.\n");   return 0; }
		
		m->mothurOut("Sampling " + toString(size) + " from each group.\n");
        
        while (lookup != NULL) {
            
            if (m->getControl_pressed()) { delete lookup; break; }
            
            processShared(lookup); delete lookup;
            
            lookup = util.getNextShared(input, allLines, userLabels, processedLabels, lastLabel);
        }
		
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getSubSampleShared");
		exit(1);
	}
}
//**********************************************************************************************************************
int SubSampleCommand::processShared(SharedRAbundVectors*& thislookup) {
	try {
        string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(sharedfile);  }
        map<string, string> variables; 
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(sharedfile));
        variables["[extension]"] = util.getExtension(sharedfile);
        variables["[distance]"] = thislookup->getLabel();
		string outputFileName = getOutputFileName("shared", variables);        
        
        SubSample sample;
        if (withReplacement)    { sample.getSampleWithReplacement(thislookup, size);  }
        else                    { sample.getSample(thislookup, size);  }
        
        if (m->getControl_pressed()) {  return 0; }
        
        ofstream out;
		util.openOutputFile(outputFileName, out);
		outputTypes["shared"].push_back(outputFileName);  outputNames.push_back(outputFileName);
        bool printHeaders = true;
		thislookup->print(out, printHeaders);
        out.close();
        
        if (constaxonomyfile != "") { //select otus from constaxonomy that are in new shared file. Also adjust the size column in the constaxonomy file
            string thisOutputDir = outputDir;
            if (outputDir == "") {  thisOutputDir += util.hasPath(constaxonomyfile);  }
            map<string, string> variables;
            variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(constaxonomyfile));
            variables["[extension]"] = util.getExtension(constaxonomyfile);
            variables["[distance]"] = thislookup->getLabel();
            string consOutputFileName = getOutputFileName("constaxonomy", variables);
            
            ifstream in;
            util.openInputFile(constaxonomyfile, in);
            
            //read headers
            string headers = util.getline(in); util.gobble(in);
            
            ofstream outCons;
            util.openOutputFile(consOutputFileName, outCons);
            outputTypes["constaxonomy"].push_back(consOutputFileName);  outputNames.push_back(consOutputFileName);

            outCons << headers << endl;
            
            while (!in.eof()) {
                
                if (m->getControl_pressed()) { break; }
                
                string otu = ""; string tax = "unknown";
                int size = 0;
                
                in >> otu; util.gobble(in); in >> size; util.gobble(in);
                tax = util.getline(in); util.gobble(in);
                
                if (m->getDebug()) { m->mothurOut("[DEBUG]: " + otu + toString(size) + tax + "\n"); }
                
                size = thislookup->getOTUTotal(otu);
                
                if (size != 0) {  outCons << otu << '\t' << size << '\t' << tax << endl; }
            }
            in.close();
            outCons.close();
        }
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "processShared");
		exit(1);
	}
}			
//**********************************************************************************************************************
int SubSampleCommand::getSubSampleList() {
	try {
        
		if (namefile != "") { util.readNames(namefile, nameMap); }
        if ((countfile == "") && withReplacement) { m->mothurOut("[WARNING]: To use the withreplacement option with a fasta file, you need to provide a count file, running without replacement.\n"); withReplacement = false;  }
        
		InputData input(listfile, "list", nullVector);
		set<string> processedLabels;
        set<string> userLabels = labels;
        string lastLabel = "";
        
        ListVector* list = util.getNextList(input, allLines, userLabels, processedLabels, lastLabel);
        
        if (constaxonomyfile != "") { //you have a constaxonomy file and have not set the labels parameter. Either your list file has one label or we only want to use one since the constaxonomy file is related to one label
            if (labels.size() == 0) {
                labels.insert(lastLabel); allLines = false;
                m->mothurOut("\n[WARNING]: The constaxonomy file represents a single label in your list file. You did not set the label parameter, so mothur is assuming you want to use label " + lastLabel + ". If this is not correct file mismatches can occur.\n\n");
            }
        }

		GroupMap groupMap;
		if (groupfile != "") {
			groupMap.readMap(groupfile);
			
			//takes care of user setting groupNames that are invalid or setting groups=all
            if (Groups.size() == 0) { Groups = groupMap.getNamesOfGroups(); }
			
			//file mismatch quit
			if (list->getNumSeqs() != groupMap.getNumSeqs()) { 
				m->mothurOut("[ERROR]: your list file contains " + toString(list->getNumSeqs()) + " sequences, and your groupfile contains " + toString(groupMap.getNumSeqs()) + ", please correct.\n"); delete list;  return 0;
			}			
		}else if (countfile != "") {
            if (ct.hasGroupInfo()) { if (Groups.size() == 0) { Groups = ct.getNamesOfGroups(); } }
            
            //file mismatch quit
			if (list->getNumSeqs() != ct.getNumUniqueSeqs()) { 
				m->mothurOut("[ERROR]: your list file contains " + toString(list->getNumSeqs()) + " sequences, and your count file contains " + toString(ct.getNumUniqueSeqs()) + " unique sequences, please correct.\n");  return 0;
			}	
        }

		//make sure that if your picked groups size is not too big
		if (persample) {
			if (size == 0) { //user has not set size, set size = smallest samples size
				if (countfile == "") { size = groupMap.getNumSeqsSmallestGroup(); }
                else {  size = ct.getNumSeqsSmallestGroup();  }
			}else { //make sure size is not too large
				vector<string> newGroups;
				for (int i = 0; i < Groups.size(); i++) {
					int thisSize = 0;
                    if (countfile == "") { thisSize = groupMap.getNumSeqs(Groups[i]); }
                    else {  thisSize = ct.getGroupCount(Groups[i]);  }
					
					if (thisSize >= size) {	newGroups.push_back(Groups[i]);	}
					else {  m->mothurOut("You have selected a size that is larger than " + Groups[i] + " number of sequences, removing " + Groups[i] + ".\n");  }
				}
				Groups = newGroups;
                if (newGroups.size() == 0) {  m->mothurOut("[ERROR]: all groups removed.\n"); m->setControl_pressed(true); }
			}
			
			m->mothurOut("Sampling " + toString(size) + " from each group.\n");
		}else{
            if (pickedGroups) {
				int total = 0;
				for(int i = 0; i < Groups.size(); i++) {
                    if (countfile == "") { total += groupMap.getNumSeqs(Groups[i]); }
                    else {  total += ct.getGroupCount(Groups[i]);  }
				}
				
				if (size == 0) { //user has not set size, set size = 10% samples size
					size = int (total * 0.10);
				}
				
				if (total < size) { 
					if (size != 0) { 
						m->mothurOut("Your size is too large for the number of groups you selected. Adjusting to " + toString(int (total * 0.10)) + ".\n");
					}
					size = int (total * 0.10);
				}
				
				m->mothurOut("Sampling " + toString(size) + " from " + toString(total) + ".\n");
			}else {
                if (size == 0) { //user has not set size, set size = 10% samples size
					if (countfile == "") {  size = int (list->getNumSeqs() * 0.10);  }
                    else { size = int (ct.getNumSeqs() * 0.10);  }
				}
				
				int thisSize = 0;
                if (countfile == "") { thisSize = list->getNumSeqs();  }
                else { thisSize = ct.getNumSeqs(); }
                
				if (size > thisSize) { m->mothurOut("Your list file only contains " + toString(thisSize) + " sequences. Setting size to " + toString(thisSize) + ".\n");
					size = thisSize;
				}
				
				m->mothurOut("Sampling " + toString(size) + " from " + toString(thisSize) + ".\n");
            }
        }
		
        set<string> subset; //dont want repeat sequence names added
		if (countfile == "") {
            SubSample sample;
            if (groupfile != "") { //use group file names
                GroupMap sampledGm = sample.getSample(groupMap, size, Groups, persample);
                vector<string> sampledSeqs = sampledGm.getNamesSeqs();
                for (int i = 0; i < sampledSeqs.size(); i++) { subset.insert(sampledSeqs[i]); }
                
                string groupOutputDir = outputDir;
                if (outputDir == "") {  groupOutputDir += util.hasPath(groupfile);  }
                map<string, string> variables;
                variables["[filename]"] = groupOutputDir + util.getRootName(util.getSimpleName(groupfile));
                variables["[extension]"] = util.getExtension(groupfile);
                string groupOutputFileName = getOutputFileName("group", variables);
                outputTypes["group"].push_back(groupOutputFileName);  outputNames.push_back(groupOutputFileName);
                sampledGm.print(groupOutputFileName);
            }else {
                vector<string> names; //use list file names
                
                //fill names
                for (int i = 0; i < list->getNumBins(); i++) {
                    string binnames = list->get(i);
                    vector<string> thisBin;
                    util.splitAtComma(binnames, thisBin);
                    
                    for(int j=0;j<thisBin.size();j++){
                        if (groupfile != "") { //if there is a groupfile given fill in group info
                            string group = groupMap.getGroup(thisBin[j]);
                            if (group == "not found") { m->mothurOut("[ERROR]: " + thisBin[j] + " is not in your groupfile. please correct.\n");  group = "NOTFOUND"; }
                            
                            //if hte user picked groups, we only want to keep the names of sequences from those groups
                            if (pickedGroups) { if (util.inUsersGroups(group, Groups)) { names.push_back(thisBin[j]); }  }
                            else{ names.push_back(thisBin[j]); }
                        }//save everyone, group
                        else{ names.push_back(thisBin[j]); }
                    }
                }
                
                util.mothurRandomShuffle(names);

                for (int i = 0; i < size; i++) { subset.insert(names[i]); }
            }
        }else {
            SubSample sample;
            CountTable sampledCt;
            if (withReplacement)    { sampledCt = sample.getSampleWithReplacement(ct, size, Groups, persample);     }
            else                    { sampledCt = sample.getSample(ct, size, Groups, persample);                    }
        
            vector<string> sampledSeqs = sampledCt.getNamesOfSeqs();
            for (int i = 0; i < sampledSeqs.size(); i++) { subset.insert(sampledSeqs[i]); }
        
            string countOutputDir = outputDir;
            if (outputDir == "") {  countOutputDir += util.hasPath(countfile);  }
            map<string, string> variables; 
            variables["[filename]"] = countOutputDir + util.getRootName(util.getSimpleName(countfile));
            variables["[extension]"] = util.getExtension(countfile);
            string countOutputFileName = getOutputFileName("count", variables);
            outputTypes["count"].push_back(countOutputFileName);  outputNames.push_back(countOutputFileName);
            sampledCt.printTable(countOutputFileName);
        }
			
        while (list != NULL) {
            
            if (m->getControl_pressed()) { delete list; break; }
            
            processList(list,  subset); delete list;
            
            list = util.getNextList(input, allLines, userLabels, processedLabels, lastLabel);
        }
		
		if (list != NULL) { delete list; }
        
        if (taxonomyfile != "") {
            if (namefile == "") {
                InputData input(listfile, "list", Groups);
                ListVector* list = input.getListVector();
                string lastLabel = list->getLabel();
                
                for (int i = 0; i < list->getNumBins(); i++) {
                    vector<string> temp;
                    string bin = list->get(i);
                    util.splitAtComma(bin, temp);
                    for (int j = 0; j < temp.size(); j++) { vector<string> tempFakeOut; tempFakeOut.push_back(temp[j]); nameMap[temp[j]] = tempFakeOut; }
                }
                delete list;
                    
                int tcount = getTax(subset);
                if (tcount != subset.size()) { m->mothurOut("[ERROR]: subsampled list file contains " + toString(subset.size()) + " sequences, but I only found " + toString(tcount) + " in your taxonomy file, did you forget a name file? Please correct.\n");  }
            }else {
                string tempAccnos = "temp.accnos";
                ofstream outAccnos;
                util.openOutputFile(tempAccnos, outAccnos);
                for (set<string>::iterator it = subset.begin(); it != subset.end(); it++) { outAccnos << *it << endl; }
                outAccnos.close();
                
                m->mothurOut("Sampling taxonomy and name file... \n");
                string thisNameOutputDir = outputDir;
                if (outputDir == "") {  thisNameOutputDir += util.hasPath(namefile);  }
                map<string, string> variables;
                variables["[filename]"] = thisNameOutputDir + util.getRootName(util.getSimpleName(namefile));
                variables["[extension]"] = util.getExtension(namefile);
                string outputNameFileName = getOutputFileName("name", variables);
                
                string thisTaxOutputDir = outputDir;
                if (outputDir == "") {  thisTaxOutputDir += util.hasPath(taxonomyfile);  }
                variables["[filename]"] = thisTaxOutputDir + util.getRootName(util.getSimpleName(taxonomyfile));
                variables["[extension]"] = util.getExtension(taxonomyfile);
                string outputTaxFileName = getOutputFileName("taxonomy", variables);
                
                
                //use unique.seqs to create new name and fastafile
                string inputString = "dups=f, name=" + namefile + ", taxonomy=" + taxonomyfile + ", accnos=" + tempAccnos;
                m->mothurOut("/******************************************/\n");
                m->mothurOut("Running command: get.seqs(" + inputString + ")\n"); 
                current->setMothurCalling(true);
                
                Command* getCommand = new GetSeqsCommand(inputString);
                getCommand->execute();
                
                map<string, vector<string> > filenames = getCommand->getOutputFiles();
                
                delete getCommand;
                current->setMothurCalling(false);
                
                util.renameFile(filenames["name"][0], outputNameFileName);
                util.renameFile(filenames["taxonomy"][0], outputTaxFileName);
                
                outputTypes["name"].push_back(outputNameFileName);  outputNames.push_back(outputNameFileName);
                outputNames.push_back(outputTaxFileName); outputTypes["taxonomy"].push_back(outputTaxFileName);
                
                m->mothurOut("/******************************************/\nDone.\n");
            }
        }
						
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getSubSampleList");
		exit(1);
	}
}
//**********************************************************************************************************************
int SubSampleCommand::processList(ListVector*& list, set<string>& subset) {
	try {
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(listfile);  }
		map<string, string> variables;
        variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(listfile));
        variables["[extension]"] = util.getExtension(listfile);
        variables["[distance]"] = list->getLabel();
		string outputFileName = getOutputFileName("list", variables);
		ofstream out;
		util.openOutputFile(outputFileName, out);
		outputTypes["list"].push_back(outputFileName);  outputNames.push_back(outputFileName);
		
		int numBins = list->getNumBins();

		ListVector* temp = new ListVector();
		temp->setLabel(list->getLabel());
		
        vector<string> binLabels = list->getLabels();
        vector<string> newLabels;
		for (int i = 0; i < numBins; i++) {
			
			if (m->getControl_pressed()) { break; }
			
			string bin = list->get(i);
            vector<string> binnames;
            util.splitAtComma(bin, binnames);
			
            string newNames = "";
			for(int j=0;j<binnames.size();j++){ if (subset.count(binnames[j]) != 0) {  newNames += binnames[j] + ",";  } }
			
			//if there are names in this bin add to new list
			if (newNames != "") { 
				newNames = newNames.substr(0, newNames.length()-1); //rip off extra comma
				temp->push_back(newNames);
                newLabels.push_back(binLabels[i]);
			}
		}
		
        temp->setLabels(newLabels);
		delete list;
		list = temp;
		
		if (m->getControl_pressed()) { out.close(); return 0; }
		
		list->print(out, false);
        out.close();
        
        if (constaxonomyfile != "") { //select otus from constaxonomy that are in new list file. Also adjust the size column in the constaxonomy file
            string thisOutputDir = outputDir;
            if (outputDir == "") {  thisOutputDir += util.hasPath(constaxonomyfile);  }
            map<string, string> variables;
            variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(constaxonomyfile));
            variables["[extension]"] = util.getExtension(constaxonomyfile);
            variables["[distance]"] = list->getLabel();
            string consOutputFileName = getOutputFileName("constaxonomy", variables);
            
            ifstream in;
            util.openInputFile(constaxonomyfile, in);
            
            //read headers
            string headers = util.getline(in); util.gobble(in);
            
            ofstream outCons;
            util.openOutputFile(consOutputFileName, outCons);
            outputTypes["constaxonomy"].push_back(consOutputFileName);  outputNames.push_back(consOutputFileName);
            
            outCons << headers << endl;
            
            while (!in.eof()) {
                
                if (m->getControl_pressed()) { break; }
                
                string otu = ""; string tax = "unknown";
                int size = 0;
                
                in >> otu; util.gobble(in); in >> size; util.gobble(in);
                tax = util.getline(in); util.gobble(in);
                
                if (m->getDebug()) { m->mothurOut("[DEBUG]: " + otu + toString(size) + tax + "\n"); }
                
                size = list->getOTUTotal(otu);
                
                if (size != 0) {  outCons << otu << '\t' << size << '\t' << tax << endl; }
            }
            in.close();
            outCons.close();
        }

		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "processList");
		exit(1);
	}
}
//**********************************************************************************************************************
int SubSampleCommand::getSubSampleRabund() {
	try {
		InputData input(rabundfile, "rabund", nullVector);
		RAbundVector* rabund = input.getRAbundVector();
		string lastLabel = rabund->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		if (size == 0) { //user has not set size, set size = 10%
			size = int((rabund->getNumSeqs()) * 0.10);
		}else if (size > rabund->getNumSeqs()) { m->mothurOut("The size you selected is too large, skipping rabund file.\n"); delete rabund; return 0; }
		
		m->mothurOut("Sampling " + toString(size) + " from " + toString(rabund->getNumSeqs()) + ".\n");
		
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(rabundfile);  }
        map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(rabundfile));
        variables["[extension]"] = util.getExtension(rabundfile);
		string outputFileName = getOutputFileName("rabund", variables);		
		ofstream out;
		util.openOutputFile(outputFileName, out);
		outputTypes["rabund"].push_back(outputFileName);  outputNames.push_back(outputFileName);
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((rabund != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			if (m->getControl_pressed()) {  delete rabund; out.close(); return 0;  }
			
			if(allLines == 1 || labels.count(rabund->getLabel()) == 1){			
				
				m->mothurOut(rabund->getLabel()+"\n");
				
				processRabund(rabund, out);
				
				processedLabels.insert(rabund->getLabel());
				userLabels.erase(rabund->getLabel());
			}
			
			if ((util.anyLabelsToProcess(rabund->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = rabund->getLabel();
				
				delete rabund; 
				
				rabund = input.getRAbundVector(lastLabel);
				m->mothurOut(rabund->getLabel()+"\n");
				
				processRabund(rabund, out);
				
				processedLabels.insert(rabund->getLabel());
				userLabels.erase(rabund->getLabel());
				
				//restore real lastlabel to save below
				rabund->setLabel(saveLabel);
			}
			
			lastLabel = rabund->getLabel();
			
			//prevent memory leak
			delete rabund; rabund = NULL;
			
			//get next line to process
			rabund = input.getRAbundVector();
		}
		
		
		if (m->getControl_pressed()) {  out.close(); return 0;  }
		
		//output error messages about any remaining user labels
		bool needToRun = false;
		for (set<string>::iterator it = userLabels.begin(); it != userLabels.end(); it++) {
			m->mothurOut("Your file does not include the label " + *it); 
            if (processedLabels.count(lastLabel) != 1)  { m->mothurOut(". I will use " + lastLabel + ".\n"); needToRun = true;  }
			else                                        { m->mothurOut(". Please refer to " + lastLabel + ".\n");               }
		}
		
		//run last label if you need to
		if (needToRun )  {
			if (rabund != NULL) { delete rabund; }
			
			rabund = input.getRAbundVector(lastLabel);
			
			m->mothurOut(rabund->getLabel()+"\n");
			
			processRabund(rabund, out);
			
			delete rabund;
		}
		out.close();  
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getSubSampleRabund");
		exit(1);
	}
}
//**********************************************************************************************************************
int SubSampleCommand::processRabund(RAbundVector*& rabund, ofstream& out) {
	try {
		
        RAbundVector* copy = new RAbundVector(*rabund);
        SubSample sample;
        if (withReplacement) {
            sample.getSampleWithReplacement(copy, size);
        }else {
            sample.getSample(copy, size);
        }
        copy->print(out);
        delete copy;
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "processRabund");
		exit(1);
	}
}	
//**********************************************************************************************************************
int SubSampleCommand::getSubSampleSabund() {
	try {
				
		InputData input(sabundfile, "sabund", nullVector);
		SAbundVector* sabund = input.getSAbundVector();
		string lastLabel = sabund->getLabel();
		
		//if the users enters label "0.06" and there is no "0.06" in their file use the next lowest label.
		set<string> processedLabels;
		set<string> userLabels = labels;
		
		if (size == 0) { //user has not set size, set size = 10%
			size = int((sabund->getNumSeqs()) * 0.10);
		}else if (size > sabund->getNumSeqs()) { m->mothurOut("The size you selected is too large, skipping sabund file.\n"); delete sabund; return 0; }
		
		
		m->mothurOut("Sampling " + toString(size) + " from " + toString(sabund->getNumSeqs()) + ".\n");
		
		string thisOutputDir = outputDir;
		if (outputDir == "") {  thisOutputDir += util.hasPath(sabundfile);  }
        map<string, string> variables; 
		variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(sabundfile));
        variables["[extension]"] = util.getExtension(sabundfile);
		string outputFileName = getOutputFileName("sabund", variables);		
		ofstream out;
		util.openOutputFile(outputFileName, out);
		outputTypes["sabund"].push_back(outputFileName);  outputNames.push_back(outputFileName);
		
		
		//as long as you are not at the end of the file or done wih the lines you want
		while((sabund != NULL) && ((allLines == 1) || (userLabels.size() != 0))) {
			if (m->getControl_pressed()) {  delete sabund; out.close(); return 0;  }
			
			if(allLines == 1 || labels.count(sabund->getLabel()) == 1){			
				
				m->mothurOut(sabund->getLabel()+"\n");
				
				processSabund(sabund, out);
				
				processedLabels.insert(sabund->getLabel());
				userLabels.erase(sabund->getLabel());
			}
			
			if ((util.anyLabelsToProcess(sabund->getLabel(), userLabels, "") ) && (processedLabels.count(lastLabel) != 1)) {
				string saveLabel = sabund->getLabel();
				
				delete sabund; 
				
				sabund = input.getSAbundVector(lastLabel);
				m->mothurOut(sabund->getLabel()+"\n");
				
				processSabund(sabund, out);
				
				processedLabels.insert(sabund->getLabel());
				userLabels.erase(sabund->getLabel());
				
				//restore real lastlabel to save below
				sabund->setLabel(saveLabel);
			}
			
			lastLabel = sabund->getLabel();
			
			//prevent memory leak
			delete sabund; sabund = NULL;
			
			//get next line to process
			sabund = input.getSAbundVector();
		}
		
		
		if (m->getControl_pressed()) {  out.close(); return 0;  }
		
		//output error messages about any remaining user labels
		bool needToRun = false;
		for (set<string>::iterator it = userLabels.begin(); it != userLabels.end(); it++) {
			m->mothurOut("Your file does not include the label " + *it); 
            if (processedLabels.count(lastLabel) != 1)  { m->mothurOut(". I will use " + lastLabel + ".\n"); needToRun = true;  }
			else                                        { m->mothurOut(". Please refer to " + lastLabel + ".\n");               }
		}
		
		//run last label if you need to
		if (needToRun )  {
			if (sabund != NULL) { delete sabund; }
			
			sabund = input.getSAbundVector(lastLabel);
			
			m->mothurOut(sabund->getLabel()+"\n");
			
			processSabund(sabund, out);
			
			delete sabund;
		}
		out.close();  
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "getSubSampleSabund");
		exit(1);
	}
}
//**********************************************************************************************************************
int SubSampleCommand::processSabund(SAbundVector*& sabund, ofstream& out) {
	try {
		SAbundVector* copy = new SAbundVector(*sabund);
        SubSample sample;
        if (withReplacement) {
            sample.getSampleWithReplacement(copy, size);
        }else {
            sample.getSample(copy, size);
        }
        copy->print(out);
        delete copy;
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SubSampleCommand", "processSabund");
		exit(1);
	}
}

//**********************************************************************************************************************
int SubSampleCommand::getTax(set<string>& subset) {
	try {

        string thisTaxOutputDir = outputDir;
        if (outputDir == "") {  thisTaxOutputDir += util.hasPath(taxonomyfile);  }
        map<string, string> variables;
        variables["[filename]"] = thisTaxOutputDir + util.getRootName(util.getSimpleName(taxonomyfile));
        variables["[extension]"] = util.getExtension(taxonomyfile);
        string outputTaxFileName = getOutputFileName("taxonomy", variables);
        ofstream outTax;
        util.openOutputFile(outputTaxFileName, outTax);
        outputNames.push_back(outputTaxFileName); outputTypes["taxonomy"].push_back(outputTaxFileName);
        
        //read through fasta file outputting only the names on the subsample list
        ifstream inTax;
        util.openInputFile(taxonomyfile, inTax);
        
        string tname, tax;
        int tcount = 0;
        map<string, vector<string> >::iterator itNameMap;
        
        while(!inTax.eof()){
            
            if (m->getControl_pressed()) { inTax.close(); outTax.close();  return 0; }
            
            inTax >> tname; util.gobble(inTax);
            tax = util.getline(inTax); util.gobble(inTax);
            
            //does the subset contain a sequence that this sequence represents
            itNameMap = nameMap.find(tname);
            if (itNameMap != nameMap.end()) {
                vector<string> nameRepresents = itNameMap->second;
                
            
                for (int i = 0; i < nameRepresents.size(); i++){
                    if (subset.count(nameRepresents[i]) != 0) {
                        outTax << nameRepresents[i] << '\t' << tax << endl;
                        tcount++;
                       
                    }
                }
            }else{ m->mothurOut("[ERROR]: " + tname + " is missing, please correct.\n");  }
        }
        inTax.close();
        outTax.close();
        
        return tcount;
    }
    catch(exception& e) {
        m->errorOut(e, "SubSampleCommand", "getTax");
        exit(1);
    }
}

//**********************************************************************************************************************



