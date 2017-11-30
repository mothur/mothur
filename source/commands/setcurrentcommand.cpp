/*
 *  setcurrentcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 3/16/11.
 *  Copyright 2011 Schloss Lab. All rights reserved.
 *
 */

#include "setcurrentcommand.h"

//**********************************************************************************************************************
vector<string> SetCurrentCommand::setParameters(){	
	try {
		
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter pflow("flow", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pflow);
        CommandParameter pfile("file", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pfile);
        CommandParameter pbiom("biom", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pbiom);
		CommandParameter pphylip("phylip", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pphylip);
		CommandParameter pcolumn("column", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pcolumn);
        CommandParameter psummary("summary", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(psummary);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pname);
		CommandParameter pgroup("group", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pgroup);
		CommandParameter plist("list", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(plist);
		CommandParameter ptaxonomy("taxonomy", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(ptaxonomy);
        CommandParameter pconstaxonomy("constaxonomy", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pconstaxonomy);
        CommandParameter pcontigsreport("contigsreport", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pcontigsreport);
		CommandParameter pqfile("qfile", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pqfile);
		CommandParameter paccnos("accnos", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(paccnos);		
		CommandParameter prabund("rabund", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(prabund);
		CommandParameter psabund("sabund", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(psabund);
		CommandParameter pdesign("design", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pdesign);
		CommandParameter porder("order", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(porder);
		CommandParameter ptree("tree", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(ptree);
		CommandParameter pshared("shared", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pshared);
		CommandParameter pordergroup("ordergroup", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pordergroup);
        CommandParameter pcount("count", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pcount);
        CommandParameter pcurrent("current", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pcurrent);
		CommandParameter prelabund("relabund", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(prelabund);
		CommandParameter psff("sff", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(psff);
		CommandParameter poligos("oligos", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(poligos);
		CommandParameter pclear("clear", "String", "", "", "", "", "","",false,false); parameters.push_back(pclear);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SetCurrentCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SetCurrentCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The set.current command allows you to set the current files saved by mothur.\n";
		helpString += "The set.current command parameters are: current, clear, phylip, column, list, rabund, sabund, name, group, design, order, tree, shared, ordergroup, relabund, fasta, qfile, sff, oligos, accnos, biom, count, summary, file, contigsreport, constaxonomy and taxonomy.\n";
        helpString += "The current parameter is used to input the output file from get.current.  This function is intended to allow you to input filenames from previous instances on mothur.  NOTE: If you have a current file set in the file *.current_files.summary file, and also set a value for that file type, the value set takes precedence.  For example, if you run set.current(current=current_files.summary, fasta=abrecovery.fasta) and your have fasta=final.fasta in the *.current_files.summary file the current fasta file will be set to abrecovery.fasta.\n";
		helpString += "The clear parameter is used to indicate which file types you would like to clear values for, multiple types can be separated by dashes.\n";
		helpString += "The set.current command should be in the following format: \n";
		helpString += "set.current(fasta=yourFastaFile) or set.current(fasta=amazon.fasta, clear=name-accnos)\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SetCurrentCommand", "getHelpString");
		exit(1);
	}
}


//**********************************************************************************************************************
SetCurrentCommand::SetCurrentCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
        vector<string> tempOutNames;
        outputTypes["summary"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SetCurrentCommand", "SetCurrentCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
string SetCurrentCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "summary") {  pattern = "[filename],current_files.summary"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "SetCurrentCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
SetCurrentCommand::SetCurrentCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			//valid paramters for this command
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
			map<string,string>::iterator it;
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
            vector<string> tempOutNames;
            outputTypes["summary"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("phylip");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["phylip"] = inputDir + it->second;		}
				}
				
				it = parameters.find("column");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["column"] = inputDir + it->second;		}
				}
				
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("list");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["list"] = inputDir + it->second;		}
				}
				
				it = parameters.find("rabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["rabund"] = inputDir + it->second;		}
				}
				
				it = parameters.find("sabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["sabund"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
				
				it = parameters.find("design");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["design"] = inputDir + it->second;		}
				}
				
				it = parameters.find("order");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["order"] = inputDir + it->second;		}
				}
				
				it = parameters.find("tree");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["tree"] = inputDir + it->second;		}
				}
				
				it = parameters.find("shared");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["shared"] = inputDir + it->second;		}
				}
				
				it = parameters.find("ordergroup");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["ordergroup"] = inputDir + it->second;		}
				}
				
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
                
				it = parameters.find("relabund");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["relabund"] = inputDir + it->second;		}
				}
				
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("qfile");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["qfile"] = inputDir + it->second;		}
				}
				
				it = parameters.find("sff");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["sff"] = inputDir + it->second;		}
				}
				
				it = parameters.find("oligos");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["oligos"] = inputDir + it->second;		}
				}
				
				it = parameters.find("accnos");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["accnos"] = inputDir + it->second;		}
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
                
                it = parameters.find("contigsreport");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["contigsreport"] = inputDir + it->second;		}
                }
				
				it = parameters.find("flow");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["flow"] = inputDir + it->second;		}
				}
                
                it = parameters.find("biom");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["biom"] = inputDir + it->second;		}
				}
                
                it = parameters.find("summary");
				//user has given a template file
				if(it != parameters.end()){
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["summary"] = inputDir + it->second;		}
				}
                
                it = parameters.find("file");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["file"] = inputDir + it->second;		}
                }
                
                it = parameters.find("current");
                //user has given a template file
                if(it != parameters.end()){
                    path = util.hasPath(it->second);
                    //if the user has not given a path then, add inputdir. else leave path alone.
                    if (path == "") {	parameters["current"] = inputDir + it->second;		}
                }
			}
            
            currentFile = validParameter.validFile(parameters, "current");
            if (currentFile == "not open") { m->mothurOut("Ignoring: " + parameters["current"]); m->mothurOutEndLine(); currentFile = ""; }
            else if (currentFile == "not found") {  currentFile = "";  }
            if (currentFile != "") { readCurrentFiles(); } //setting variables overwrites the settings in the file.
            
            //if the user changes the output directory command factory will send this info to us in the output parameter
            outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){  outputDir = ""; }
			
			//check for parameters
			phylipfile = validParameter.validFile(parameters, "phylip");
			if (phylipfile == "not open") { m->mothurOut("Ignoring: " + parameters["phylip"]); m->mothurOutEndLine(); phylipfile = ""; }
			else if (phylipfile == "not found") {  phylipfile = "";  }	
			if (phylipfile != "") { current->setPhylipFile(phylipfile); }
			
			columnfile = validParameter.validFile(parameters, "column");
			if (columnfile == "not open") { m->mothurOut("Ignoring: " + parameters["column"]); m->mothurOutEndLine(); columnfile = ""; }
			else if (columnfile == "not found") {  columnfile = "";  }	
			if (columnfile != "") { current->setColumnFile(columnfile); }
			
			listfile = validParameter.validFile(parameters, "list");
			if (listfile == "not open") { m->mothurOut("Ignoring: " + parameters["list"]); m->mothurOutEndLine(); listfile = ""; }
			else if (listfile == "not found") {  listfile = "";  }	
			if (listfile != "") { current->setListFile(listfile); }
			
			rabundfile = validParameter.validFile(parameters, "rabund");
			if (rabundfile == "not open") { m->mothurOut("Ignoring: " + parameters["rabund"]); m->mothurOutEndLine(); rabundfile = ""; }
			else if (rabundfile == "not found") {  rabundfile = "";  }	
			if (rabundfile != "") { current->setRabundFile(rabundfile); }
			
			sabundfile = validParameter.validFile(parameters, "sabund");
			if (sabundfile == "not open") { m->mothurOut("Ignoring: " + parameters["sabund"]); m->mothurOutEndLine(); sabundfile = ""; }
			else if (sabundfile == "not found") {  sabundfile = "";  }	
			if (sabundfile != "") { current->setSabundFile(sabundfile); }
			
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { m->mothurOut("Ignoring: " + parameters["name"]); m->mothurOutEndLine(); namefile = ""; }
			else if (namefile == "not found") {  namefile = "";  }	
			if (namefile != "") { current->setNameFile(namefile); }
			
			groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not open") { m->mothurOut("Ignoring: " + parameters["group"]); m->mothurOutEndLine(); groupfile = ""; }
			else if (groupfile == "not found") {  groupfile = "";  }
			if (groupfile != "") { current->setGroupFile(groupfile); }
            
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { m->mothurOut("Ignoring: " + parameters["count"]); m->mothurOutEndLine(); countfile = ""; }
			else if (countfile == "not found") {  countfile = "";  }
			if (countfile != "") { current->setCountFile(countfile); }
			
			designfile = validParameter.validFile(parameters, "design");
			if (designfile == "not open") { m->mothurOut("Ignoring: " + parameters["design"]); m->mothurOutEndLine(); designfile = ""; }
			else if (designfile == "not found") {  designfile = "";  }	
			if (designfile != "") { current->setDesignFile(designfile); }
			
			orderfile = validParameter.validFile(parameters, "order");
			if (orderfile == "not open") { m->mothurOut("Ignoring: " + parameters["order"]); m->mothurOutEndLine(); orderfile = ""; }
			else if (orderfile == "not found") {  orderfile = "";  }
			if (orderfile != "") { current->setOrderFile(orderfile); }
			
			treefile = validParameter.validFile(parameters, "tree");
			if (treefile == "not open") { m->mothurOut("Ignoring: " + parameters["tree"]); m->mothurOutEndLine(); treefile = ""; }
			else if (treefile == "not found") {  treefile = "";  }	
			if (treefile != "") { current->setTreeFile(treefile); }
			
			sharedfile = validParameter.validFile(parameters, "shared");
			if (sharedfile == "not open") { m->mothurOut("Ignoring: " + parameters["shared"]); m->mothurOutEndLine(); sharedfile = ""; }
			else if (sharedfile == "not found") {  sharedfile = "";  }	
			if (sharedfile != "") { current->setSharedFile(sharedfile); }
			
			ordergroupfile = validParameter.validFile(parameters, "ordergroup");
			if (ordergroupfile == "not open") { m->mothurOut("Ignoring: " + parameters["ordergroup"]); m->mothurOutEndLine(); ordergroupfile = ""; }
			else if (ordergroupfile == "not found") {  ordergroupfile = "";  }	
			if (ordergroupfile != "") { current->setOrderGroupFile(ordergroupfile); }
			
			relabundfile = validParameter.validFile(parameters, "relabund");
			if (relabundfile == "not open") { m->mothurOut("Ignoring: " + parameters["relabund"]); m->mothurOutEndLine(); relabundfile = ""; }
			else if (relabundfile == "not found") {  relabundfile = "";  }	
			if (relabundfile != "") { current->setRelAbundFile(relabundfile); }
			
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not open") { m->mothurOut("Ignoring: " + parameters["fasta"]); m->mothurOutEndLine(); fastafile = ""; }
			else if (fastafile == "not found") {  fastafile = "";  }
			if (fastafile != "") { current->setFastaFile(fastafile); }
			
			qualfile = validParameter.validFile(parameters, "qfile");
			if (qualfile == "not open") { m->mothurOut("Ignoring: " + parameters["qfile"]); m->mothurOutEndLine(); qualfile = ""; }
			else if (qualfile == "not found") {  qualfile = "";  }	
			if (qualfile != "") { current->setQualFile(qualfile); }
			
			sfffile = validParameter.validFile(parameters, "sff");
			if (sfffile == "not open") { m->mothurOut("Ignoring: " + parameters["sff"]); m->mothurOutEndLine(); sfffile = ""; }
			else if (sfffile == "not found") {  sfffile = "";  }	
			if (sfffile != "") { current->setSFFFile(sfffile); }
			
			oligosfile = validParameter.validFile(parameters, "oligos");
			if (oligosfile == "not open") { m->mothurOut("Ignoring: " + parameters["oligos"]); m->mothurOutEndLine(); oligosfile = ""; }
			else if (oligosfile == "not found") {  oligosfile = "";  }	
			if (oligosfile != "") { current->setOligosFile(oligosfile); }
			
			accnosfile = validParameter.validFile(parameters, "accnos");
			if (accnosfile == "not open") { m->mothurOut("Ignoring: " + parameters["accnos"]); m->mothurOutEndLine(); accnosfile = ""; }
			else if (accnosfile == "not found") {  accnosfile = "";  }	
			if (accnosfile != "") { current->setAccnosFile(accnosfile); }
			
			taxonomyfile = validParameter.validFile(parameters, "taxonomy");
			if (taxonomyfile == "not open") { m->mothurOut("Ignoring: " + parameters["taxonomy"]); m->mothurOutEndLine(); taxonomyfile = ""; }
			else if (taxonomyfile == "not found") {  taxonomyfile = "";  }	
			if (taxonomyfile != "") { current->setTaxonomyFile(taxonomyfile); }
            
            contigsreportfile = validParameter.validFile(parameters, "contigsreport");
            if (contigsreportfile == "not open") { m->mothurOut("Ignoring: " + parameters["contigsreport"]); m->mothurOutEndLine(); contigsreportfile = ""; }
            else if (contigsreportfile == "not found") {  contigsreportfile = "";  }
            if (contigsreportfile != "") { current->setContigsReportFile(contigsreportfile); }
            
            constaxonomyfile = validParameter.validFile(parameters, "constaxonomy");
            if (constaxonomyfile == "not open") { m->mothurOut("Ignoring: " + parameters["constaxonomy"]); m->mothurOutEndLine(); constaxonomyfile = ""; }
            else if (constaxonomyfile == "not found") {  constaxonomyfile = "";  }
            if (constaxonomyfile != "") { current->setConsTaxonomyFile(constaxonomyfile); }
			
			flowfile = validParameter.validFile(parameters, "flow");
			if (flowfile == "not open") { m->mothurOut("Ignoring: " + parameters["flow"]); m->mothurOutEndLine(); flowfile = ""; }
			else if (flowfile == "not found") {  flowfile = "";  }	
			if (flowfile != "") { current->setFlowFile(flowfile); }
            
            biomfile = validParameter.validFile(parameters, "biom");
			if (biomfile == "not open") { m->mothurOut("Ignoring: " + parameters["biom"]); m->mothurOutEndLine(); biomfile = ""; }
			else if (biomfile == "not found") {  biomfile = "";  }	
			if (biomfile != "") { current->setBiomFile(biomfile); }
            
            summaryfile = validParameter.validFile(parameters, "summary");
			if (summaryfile == "not open") { m->mothurOut("Ignoring: " + parameters["summary"]); m->mothurOutEndLine(); summaryfile = ""; }
			else if (summaryfile == "not found") {  summaryfile = "";  }
			if (summaryfile != "") { current->setSummaryFile(summaryfile); }
            
            filefile = validParameter.validFile(parameters, "file");
            if (filefile == "not open") { m->mothurOut("Ignoring: " + parameters["file"]); m->mothurOutEndLine(); filefile = ""; }
            else if (filefile == "not found") {  filefile = "";  }
            if (filefile != "") { current->setFileFile(filefile); }

			string temp = validParameter.valid(parameters, "processors");
			if (temp == "not found"){	temp = current->getProcessors();	}
			current->setProcessors(temp);
			
			clearTypes = validParameter.valid(parameters, "clear");
			if (clearTypes == "not found") { clearTypes = ""; }
			else { util.splitAtDash(clearTypes, types);	}
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "SetCurrentCommand", "SetCurrentCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int SetCurrentCommand::execute(){
	try {
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
        
		//user wants to clear a type
		if (types.size() != 0) {
			for (int i = 0; i < types.size(); i++) {
				
				if (m->getControl_pressed()) { break; }
				
				//look for file types
				if (types[i] == "fasta") {
					current->setFastaFile("");
				}else if (types[i] == "qfile") {
					current->setQualFile("");
				}else if (types[i] == "phylip") {
					current->setPhylipFile("");
				}else if (types[i] == "column") {
					current->setColumnFile("");
				}else if (types[i] == "list") {
					current->setListFile("");
				}else if (types[i] == "rabund") {
					current->setRabundFile("");
				}else if (types[i] == "sabund") {
					current->setSabundFile("");
				}else if (types[i] == "name") {
					current->setNameFile("");
				}else if (types[i] == "group") {
					current->setGroupFile("");
				}else if (types[i] == "order") {
					current->setOrderFile("");
				}else if (types[i] == "ordergroup") {
					current->setOrderGroupFile("");
				}else if (types[i] == "tree") {
					current->setTreeFile("");
				}else if (types[i] == "shared") {
					current->setSharedFile("");
				}else if (types[i] == "relabund") {
					current->setRelAbundFile("");
				}else if (types[i] == "design") {
					current->setDesignFile("");
				}else if (types[i] == "sff") {
					current->setSFFFile("");
				}else if (types[i] == "oligos") {
					current->setOligosFile("");
				}else if (types[i] == "accnos") {
					current->setAccnosFile("");
				}else if (types[i] == "taxonomy") {
					current->setTaxonomyFile("");
                }else if (types[i] == "constaxonomy") {
                    current->setConsTaxonomyFile("");
                }else if (types[i] == "contigsreport") {
                    current->setContigsReportFile("");
				}else if (types[i] == "flow") {
					current->setFlowFile("");
                }else if (types[i] == "biom") {
					current->setBiomFile("");
                }else if (types[i] == "count") {
					current->setCountFile("");
                }else if (types[i] == "summary") {
					current->setSummaryFile("");
                }else if (types[i] == "file") {
                    current->setFileFile("");
				}else if (types[i] == "processors") {
					current->setProcessors("1");
				}else if (types[i] == "all") {
					current->clearCurrentFiles();
				}else {
					m->mothurOut("[ERROR]: mothur does not save a current file for " + types[i] + "\n");
				}
			}
		}
        
		m->mothurOutEndLine(); m->mothurOut("Current files saved by mothur:"); m->mothurOutEndLine();
        
        if (current->hasCurrentFiles()) {
            map<string, string> variables;
            variables["[filename]"] = util.getFullPathName(outputDir);
            string filename = getOutputFileName("summary", variables);
            
            current->printCurrentFiles(filename);
            outputNames.push_back(filename); outputTypes["summary"].push_back(filename);
            
            m->mothurOutEndLine();
            m->mothurOut("Output File Names: "); m->mothurOutEndLine();
            for (int i = 0; i < outputNames.size(); i++) { m->mothurOut(outputNames[i]); m->mothurOutEndLine(); }
            m->mothurOutEndLine();
        }
		
		return 0;	
	}
	
	catch(exception& e) {
		m->errorOut(e, "SetCurrentCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int SetCurrentCommand::readCurrentFiles(){
    try{
        
        ifstream in;
        util.openInputFile(currentFile, in);
        
        
        while(!in.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            string line = util.getline(in); util.gobble(in);
            
            vector<string> pieces;
            util.splitAtChar(line, pieces, '=');
            
            if (pieces.size() != 2) { m->mothurOut("[ERROR]: " + util.getStringFromVector(pieces, ",") + " line is not in the correct format.  Did you edit the file? Mothur expects tag=filename.  Example: fasta=final.fasta\n"); m->setControl_pressed(true);  }
            else{
                   //look for file types
                if (pieces[0] == "fasta") {
                    current->setFastaFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "qfile") {
                    current->setQualFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "phylip") {
                    current->setPhylipFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "column") {
                    current->setColumnFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "list") {
                    current->setListFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "rabund") {
                    current->setRabundFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "sabund") {
                    current->setSabundFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "name") {
                    current->setNameFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "group") {
                    current->setGroupFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "order") {
                    current->setOrderFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "ordergroup") {
                    current->setOrderGroupFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "tree") {
                    current->setTreeFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "shared") {
                    current->setSharedFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "relabund") {
                    current->setRelAbundFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "design") {
                    current->setDesignFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "sff") {
                    current->setSFFFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "oligos") {
                    current->setOligosFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "accnos") {
                    current->setAccnosFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "taxonomy") {
                    current->setTaxonomyFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "flow") {
                    current->setFlowFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "biom") {
                    current->setBiomFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "count") {
                    current->setCountFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "summary") {
                    current->setSummaryFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "file") {
                    current->setFileFile(util.getFullPathName(pieces[1]));
                }else if (pieces[0] == "processors") {
                    current->setProcessors(pieces[1]);
                }else {
                    m->mothurOut("[ERROR]: mothur does not save a current file for " + util.getFullPathName(pieces[1])); m->mothurOutEndLine();
                }
            }
        }
        in.close();
        
        return 0;
    }
    
    catch(exception& e) {
        m->errorOut(e, "SetCurrentCommand", "readCurrentFiles");
        exit(1);
    }
}
//**********************************************************************************************************************



