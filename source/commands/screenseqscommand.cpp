/*
 *  screenseqscommand.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 6/3/09.
 *  Copyright 2009 Patrick D. Schloss. All rights reserved.
 *
 */

#include "screenseqscommand.h"
#include "counttable.h"
#include "summary.hpp"
#include "removeseqscommand.h"

//**********************************************************************************************************************
vector<string> ScreenSeqsCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","fasta",false,true,true); parameters.push_back(pfasta);
        CommandParameter pcontigsreport("contigsreport", "InputTypes", "", "", "report", "none", "none","contigsreport",false,false,true); parameters.push_back(pcontigsreport);
        CommandParameter palignreport("alignreport", "InputTypes", "", "", "report", "none", "none","alignreport",false,false); parameters.push_back(palignreport);
        CommandParameter psummary("summary", "InputTypes", "", "", "report", "none", "none","summary",false,false); parameters.push_back(psummary);
        CommandParameter pname("name", "InputTypes", "", "", "NameCount", "none", "none","name",false,false,true); parameters.push_back(pname);
        CommandParameter pcount("count", "InputTypes", "", "", "NameCount-CountGroup", "none", "none","count",false,false,true); parameters.push_back(pcount);
		CommandParameter pgroup("group", "InputTypes", "", "", "CountGroup", "none", "none","group",false,false,true); parameters.push_back(pgroup);
		CommandParameter pqfile("qfile", "InputTypes", "", "", "none", "none", "none","qfile",false,false); parameters.push_back(pqfile);
		
		CommandParameter ptax("taxonomy", "InputTypes", "", "", "none", "none", "none","taxonomy",false,false); parameters.push_back(ptax);
		CommandParameter pstart("start", "Number", "", "-1", "", "", "","",false,false,true); parameters.push_back(pstart);
		CommandParameter pend("end", "Number", "", "-1", "", "", "","",false,false,true); parameters.push_back(pend);
		CommandParameter pmaxambig("maxambig", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pmaxambig);
		CommandParameter pmaxhomop("maxhomop", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pmaxhomop);
		CommandParameter pminlength("minlength", "Number", "", "10", "", "", "","",false,false); parameters.push_back(pminlength);
		CommandParameter pmaxlength("maxlength", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pmaxlength);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter pcriteria("criteria", "Number", "", "90", "", "", "","",false,false); parameters.push_back(pcriteria);
		CommandParameter poptimize("optimize", "Multiple", "none-start-end-maxambig-maxhomop-minlength-maxlength", "none", "", "", "","",true,false); parameters.push_back(poptimize);
        CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        //report parameters
        CommandParameter pminoverlap("minoverlap", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pminoverlap);
        CommandParameter postart("ostart", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(postart);
        CommandParameter poend("oend", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(poend);
        CommandParameter pmismatches("mismatches", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pmismatches);
        CommandParameter pmaxn("maxn", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pmaxn);
        CommandParameter pminscore("minscore", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pminscore);
        CommandParameter pmaxinsert("maxinsert", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pmaxinsert);
        CommandParameter pminsim("minsim", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pminsim);

		
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ScreenSeqsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The screen.seqs command reads a fastafile and screens sequences.\n";
		helpString += "The screen.seqs command parameters are fasta, start, end, maxambig, maxhomop, minlength, maxlength, name, group, count, qfile, alignreport, contigsreport, summary, taxonomy, optimize, criteria and processors.\n";
		helpString += "The fasta parameter is required.\n";
        helpString += "The contigsreport parameter allows you to use the contigsreport file to determine if a sequence is good. Screening parameters include: minoverlap, ostart, oend and mismatches. \n";
        helpString += "The alignreport parameter allows you to use the alignreport file to determine if a sequence is good. Screening parameters include: minsim, minscore and maxinsert. \n";
        helpString += "The summary parameter allows you to use the summary file from summary.seqs to save time processing.\n";
		helpString += "The taxonomy parameter allows you to remove bad seqs from taxonomy files.\n";
		helpString += "The start parameter is used to set a position the \"good\" sequences must start by. The default is -1.\n";
		helpString += "The end parameter is used to set a position the \"good\" sequences must end after. The default is -1.\n";
		helpString += "The maxambig parameter allows you to set the maximum number of ambiguous bases allowed. The default is -1.\n";
		helpString += "The maxhomop parameter allows you to set a maximum homopolymer length. \n";
		helpString += "The minlength parameter allows you to set and minimum sequence length. Default=10.\n";
		helpString += "The maxn parameter allows you to set and maximum number of N's allowed in a sequence. \n";
        helpString += "The minoverlap parameter allows you to set and minimum overlap. The default is -1. \n";
        helpString += "The ostart parameter is used to set an overlap position the \"good\" sequences must start by. The default is -1. \n";
        helpString += "The oend parameter is used to set an overlap position the \"good\" sequences must end after. The default is -1.\n";
        helpString += "The mismatches parameter allows you to set and maximum mismatches in the contigs.report. \n";
        helpString += "The minsim parameter allows you to set the minimum similarity to template sequences during alignment. Found in column \'SimBtwnQuery&Template\' in align.report file.\n";
        helpString += "The minscore parameter allows you to set the minimum search score during alignment. Found in column \'SearchScore\' in align.report file.\n";
        helpString += "The maxinsert parameter allows you to set the maximum number of insertions during alignment. Found in column \'LongestInsert\' in align.report file.\n";
		helpString += "The processors parameter allows you to specify the number of processors to use while running the command. The default is 1.\n";
		helpString += "The optimize and criteria parameters allow you set the start, end, maxabig, maxhomop, minlength and maxlength parameters relative to your set of sequences .\n";
		helpString += "For example optimize=start-end, criteria=90, would set the start and end values to the position 90% of your sequences started and ended.\n";
		helpString += "The name parameter allows you to provide a namesfile, and the group parameter allows you to provide a groupfile.\n";
		helpString += "The screen.seqs command should be in the following format: \n";
		helpString += "screen.seqs(fasta=yourFastaFile, name=youNameFile, group=yourGroupFIle, start=yourStart, end=yourEnd, maxambig=yourMaxambig,  \n";
		helpString += "maxhomop=yourMaxhomop, minlength=youMinlength, maxlength=yourMaxlength)  \n";	
		helpString += "Example screen.seqs(fasta=abrecovery.fasta, name=abrecovery.names, group=abrecovery.groups, start=..., end=..., maxambig=..., maxhomop=..., minlength=..., maxlength=...).\n";
		;
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ScreenSeqsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta")            {   pattern = "[filename],good,[extension]";    }
        else if (type == "taxonomy")    {   pattern = "[filename],good,[extension]";    }
        else if (type == "name")        {   pattern = "[filename],good,[extension]";    }
        else if (type == "group")       {   pattern = "[filename],good,[extension]";    }
        else if (type == "count")       {   pattern = "[filename],good,[extension]";    }
        else if (type == "accnos")      {   pattern = "[filename],bad.accnos";          }
        else if (type == "qfile")       {   pattern = "[filename],good,[extension]";    }
        else if (type == "alignreport")      {   pattern = "[filename],good.align.report";    }
        else if (type == "contigsreport")      {   pattern = "[filename],good.contigs.report";    }
        else if (type == "summary")      {   pattern = "[filename],good.summary";    }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ScreenSeqsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
ScreenSeqsCommand::ScreenSeqsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
		outputTypes["name"] = tempOutNames;
		outputTypes["group"] = tempOutNames;
		outputTypes["alignreport"] = tempOutNames;
        outputTypes["contigsreport"] = tempOutNames;
        outputTypes["summary"] = tempOutNames;
		outputTypes["accnos"] = tempOutNames;
		outputTypes["qfile"] = tempOutNames;
		outputTypes["taxonomy"] = tempOutNames;
        outputTypes["count"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "ScreenSeqsCommand");
		exit(1);
	}
}
//***************************************************************************************************************

ScreenSeqsCommand::ScreenSeqsCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter("screen.seqs");
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["fasta"] = tempOutNames;
			outputTypes["name"] = tempOutNames;
			outputTypes["group"] = tempOutNames;
			outputTypes["alignreport"] = tempOutNames;
			outputTypes["accnos"] = tempOutNames;
			outputTypes["qfile"] = tempOutNames;
			outputTypes["taxonomy"] = tempOutNames;
            outputTypes["count"] = tempOutNames;
			outputTypes["contigsreport"] = tempOutNames;
            outputTypes["summary"] = tempOutNames;

            
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
				
				it = parameters.find("alignreport");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["alignreport"] = inputDir + it->second;		}
				}
                
                it = parameters.find("contigsreport");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["contigsreport"] = inputDir + it->second;		}
				}
                
                it = parameters.find("summary");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["summary"] = inputDir + it->second;		}
				}
				
				it = parameters.find("qfile");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["qfile"] = inputDir + it->second;		}
				}
				
				it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
			}

            fileType = "name file";
			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta");
			if (fastafile == "not found") { 			
				fastafile = current->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (fastafile == "not open") { abort = true; }
			else { current->setFastaFile(fastafile); }
	
			groupfile = validParameter.validFile(parameters, "group");
			if (groupfile == "not open") { abort = true; }	
			else if (groupfile == "not found") { groupfile = ""; }
			else { current->setGroupFile(groupfile); }
			
			qualfile = validParameter.validFile(parameters, "qfile");
			if (qualfile == "not open") { abort = true; }	
			else if (qualfile == "not found") { qualfile = ""; }
			else { current->setQualFile(qualfile); }
			
			namefile = validParameter.validFile(parameters, "name");
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") { namefile = ""; }	
			else { current->setNameFile(namefile); }
			
            countfile = validParameter.validFile(parameters, "count");
			if (countfile == "not open") { countfile = ""; abort = true; }
			else if (countfile == "not found") { countfile = "";  }	
			else { current->setCountFile(countfile); fileType = "count file"; }
            
            contigsreport = validParameter.validFile(parameters, "contigsreport");
			if (contigsreport == "not open") { contigsreport = ""; abort = true; }
			else if (contigsreport == "not found") { contigsreport = "";  }	
            
            summaryfile = validParameter.validFile(parameters, "summary");
			if (summaryfile == "not open") { summaryfile = ""; abort = true; }
			else if (summaryfile == "not found") { summaryfile = "";  }
            else { current->setSummaryFile(summaryfile); }
            
            if ((namefile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: name or count."); m->mothurOutEndLine(); abort = true;
            }
			
            if ((groupfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: group or count."); m->mothurOutEndLine(); abort=true;
            }
            
			alignreport = validParameter.validFile(parameters, "alignreport");
			if (alignreport == "not open") { abort = true; }
			else if (alignreport == "not found") { alignreport = ""; }
			
			taxonomy = validParameter.validFile(parameters, "taxonomy");
			if (taxonomy == "not open") { abort = true; }
			else if (taxonomy == "not found") { taxonomy = ""; }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += util.hasPath(fastafile); //if user entered a file with a path then preserve it	
			}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.valid(parameters, "start");		if (temp == "not found") { temp = "-1"; }
			util.mothurConvert(temp, startPos); 
		
			temp = validParameter.valid(parameters, "end");			if (temp == "not found") { temp = "-1"; }
			util.mothurConvert(temp, endPos);  

			temp = validParameter.valid(parameters, "maxambig");		if (temp == "not found") { temp = "-1"; }
			util.mothurConvert(temp, maxAmbig);  

			temp = validParameter.valid(parameters, "maxhomop");		if (temp == "not found") { temp = "-1"; }
			util.mothurConvert(temp, maxHomoP);  

			temp = validParameter.valid(parameters, "minlength");	if (temp == "not found") { temp = "10"; }
			util.mothurConvert(temp, minLength); 
			
			temp = validParameter.valid(parameters, "maxlength");	if (temp == "not found") { temp = "-1"; }
			util.mothurConvert(temp, maxLength); 
			
			temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);
			
            temp = validParameter.valid(parameters, "minoverlap");	if (temp == "not found") { temp = "-1"; }
			util.mothurConvert(temp, minOverlap); 
            
            temp = validParameter.valid(parameters, "ostart");	if (temp == "not found") { temp = "-1"; }
			util.mothurConvert(temp, oStart); 
            
            temp = validParameter.valid(parameters, "oend");	if (temp == "not found") { temp = "-1"; }
			util.mothurConvert(temp, oEnd); 
            
            temp = validParameter.valid(parameters, "mismatches");	if (temp == "not found") { temp = "-1"; }
			util.mothurConvert(temp, mismatches); 
            
            temp = validParameter.valid(parameters, "maxn");	if (temp == "not found") { temp = "-1"; }
			util.mothurConvert(temp, maxN); 
            
            temp = validParameter.valid(parameters, "minscore");	if (temp == "not found") { temp = "-1"; }
			util.mothurConvert(temp, minScore); 
            
            temp = validParameter.valid(parameters, "maxinsert");	if (temp == "not found") { temp = "-1"; }
			util.mothurConvert(temp, maxInsert); 
            
            temp = validParameter.valid(parameters, "minsim");	if (temp == "not found") { temp = "-1"; }
			util.mothurConvert(temp, minSim); 
            
			temp = validParameter.valid(parameters, "optimize");	//optimizing trumps the optimized values original value
			if (temp == "not found"){	temp = "none";		}
			util.splitAtDash(temp, optimize);		
            
            if ((contigsreport != "") && ((summaryfile != "") || ( alignreport != ""))) {
                m->mothurOut("[ERROR]: You may only provide one of the following: contigsreport, alignreport or summary, aborting.\n"); abort=true;
            }
            
            if ((alignreport != "") && ((summaryfile != "") || ( contigsreport != ""))) {
                m->mothurOut("[ERROR]: You may only provide one of the following: contigsreport, alignreport or summary, aborting.\n"); abort=true;
            }
            
            if ((summaryfile != "") && ((alignreport != "") || ( contigsreport != ""))) {
                m->mothurOut("[ERROR]: You may only provide one of the following: contigsreport, alignreport or summary, aborting.\n"); abort=true;
            }
			
            //check to make sure you have the files you need for certain screening
            if ((contigsreport == "") && ((minOverlap != -1) || (oStart != -1) || (oEnd != -1) || (mismatches != -1))) {
                m->mothurOut("[ERROR]: minoverlap, ostart, oend and mismatches can only be used with a contigs.report file, aborting.\n"); abort=true;
            }
            
            if ((alignreport == "") && ((minScore != -1) || (maxInsert != -1) || (minSim != -1))) {
                m->mothurOut("[ERROR]: minscore, maxinsert and minsim can only be used with a align.report file, aborting.\n"); abort=true;
            }
            
			//check for invalid optimize options
			set<string> validOptimizers;
			validOptimizers.insert("none"); validOptimizers.insert("start"); validOptimizers.insert("end"); validOptimizers.insert("maxambig"); validOptimizers.insert("maxhomop"); validOptimizers.insert("minlength"); validOptimizers.insert("maxlength"); validOptimizers.insert("maxn");
            if (contigsreport != "")    { validOptimizers.insert("minoverlap"); validOptimizers.insert("ostart"); validOptimizers.insert("oend"); validOptimizers.insert("mismatches");  }
            if (alignreport != "")      { validOptimizers.insert("minscore"); validOptimizers.insert("maxinsert"); validOptimizers.insert("minsim"); }
            
			for (int i = 0; i < optimize.size(); i++) { 
				if (validOptimizers.count(optimize[i]) == 0) { 
					m->mothurOut(optimize[i] + " is not a valid optimizer with your input files. Valid options are "); 
                    string valid = "";
                    for (set<string>::iterator it = validOptimizers.begin(); it != validOptimizers.end(); it++) {
                        valid += (*it) + ", ";
                    }
                    if (valid.length() != 0) {  valid = valid.substr(0, valid.length()-2); }
                    m->mothurOut(valid + ".");
                    m->mothurOutEndLine();
					optimize.erase(optimize.begin()+i);
					i--;
				}
			}
			
			if (optimize.size() == 1) { if (optimize[0] == "none") { optimize.clear(); } }
			
			temp = validParameter.valid(parameters, "criteria");	if (temp == "not found"){	temp = "90";				}
			util.mothurConvert(temp, criteria); 
			
			if (countfile == "") { 
                if (namefile == "") {
                    vector<string> files; files.push_back(fastafile);
                    if (!current->getMothurCalling())  {  parser.getNameFile(files);  }
                }
            }
		}

	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "ScreenSeqsCommand");
		exit(1);
	}
}
//***************************************************************************************************************

int ScreenSeqsCommand::execute(){
	try{
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
        map<string, string> badSeqNames;
        long start = time(NULL);
        long long numFastaSeqs = 0;
        
        //use the namefile to optimize correctly
        if (namefile != "") { nameMap = util.readNames(namefile); }
        else if (countfile != "") {
            CountTable ct;
            ct.readTable(countfile, true, false);
            nameMap = ct.getNameMap();
        }
        
        map<string, string> variables;
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(fastafile));
        badAccnosFile =  getOutputFileName("accnos",variables);
        outputNames.push_back(badAccnosFile); outputTypes["accnos"].push_back(badAccnosFile);

        if ((contigsreport == "") && (summaryfile == "") && (alignreport == "")) {   numFastaSeqs = screenFasta(badSeqNames);  }
        else {   numFastaSeqs = screenReports(badSeqNames);   }
		
        if (m->getControl_pressed()) {  for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]); } return 0; }
        
        //use remove.seqs to create new name, group and count file
        if ((countfile != "") || (namefile != "") || (groupfile != "") || (qualfile != "") || (taxonomy != "")) {
            string inputString = "accnos=" + badAccnosFile;
            
            if (countfile != "") {  inputString += ", count=" + countfile;  }
            else{
                if (namefile != "") {  inputString += ", name=" + namefile;  }
                if (groupfile != "") {  inputString += ", group=" + groupfile;  }
            }
            if(qualfile != "")						{	inputString += ", qfile=" + qualfile;       }
            if(taxonomy != "")						{	inputString += ", taxonomy=" + taxonomy;	}
            
            m->mothurOut("/******************************************/"); m->mothurOutEndLine();
            m->mothurOut("Running command: remove.seqs(" + inputString + ")"); m->mothurOutEndLine();
            current->setMothurCalling(true);
            
            Command* removeCommand = new RemoveSeqsCommand(inputString);
            removeCommand->execute();
            
            map<string, vector<string> > filenames = removeCommand->getOutputFiles();
            
            delete removeCommand;
            current->setMothurCalling(false);
            m->mothurOut("/******************************************/"); m->mothurOutEndLine();
            
            if (groupfile != "") {
                string thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += util.hasPath(groupfile);  }
                variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(groupfile));
                variables["[extension]"] = util.getExtension(groupfile);
                string outGroup = getOutputFileName("group", variables);
                util.renameFile(filenames["group"][0], outGroup);
                outputNames.push_back(outGroup); outputTypes["group"].push_back(outGroup);
            }
            
            if (namefile != "") {
                string thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += util.hasPath(namefile);  }
                variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(namefile));
                variables["[extension]"] = util.getExtension(namefile);
                string outName = getOutputFileName("name", variables);
                util.renameFile(filenames["name"][0], outName);
                outputNames.push_back(outName); outputTypes["name"].push_back(outName);
            }
            
            if (countfile != "") {
                string thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += util.hasPath(countfile);  }
                variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(countfile));
                variables["[extension]"] = util.getExtension(countfile);
                string outCount = getOutputFileName("count", variables);
                util.renameFile(filenames["count"][0], outCount);
                outputNames.push_back(outCount); outputTypes["count"].push_back(outCount);
            }
            
            if (qualfile != "") {
                string thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += util.hasPath(qualfile);  }
                variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(qualfile));
                variables["[extension]"] = util.getExtension(qualfile);
                string outQual = getOutputFileName("qfile", variables);
                util.renameFile(filenames["qfile"][0], outQual);
                outputNames.push_back(outQual); outputTypes["name"].push_back(outQual);
            }
            
            if (taxonomy != "") {
                string thisOutputDir = outputDir;
                if (outputDir == "") {  thisOutputDir += util.hasPath(taxonomy);  }
                variables["[filename]"] = thisOutputDir + util.getRootName(util.getSimpleName(taxonomy));
                variables["[extension]"] = util.getExtension(taxonomy);
                string outTax = getOutputFileName("taxonomy", variables);
                util.renameFile(filenames["taxonomy"][0], outTax);
                outputNames.push_back(outTax); outputTypes["count"].push_back(outTax);
            }
        }
		
		if (m->getControl_pressed()) {  for (int i = 0; i < outputNames.size(); i++) { util.mothurRemove(outputNames[i]);  } return 0; }

        m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) { m->mothurOut(outputNames[i]); m->mothurOutEndLine(); }
		m->mothurOutEndLine();
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
		
		itTypes = outputTypes.find("qfile");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setQualFile(currentName); }
		}
		
		itTypes = outputTypes.find("taxonomy");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setTaxonomyFile(currentName); }
		}
        
        itTypes = outputTypes.find("count");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setCountFile(currentName); }
		}
        
        itTypes = outputTypes.find("accnos");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setAccnosFile(currentName); }
        }

		m->mothurOut("It took " + toString(time(NULL) - start) + " secs to screen " + toString(numFastaSeqs) + " sequences.");
		m->mothurOutEndLine();

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "execute");
		exit(1);
	}
}
//***************************************************************************************************************/
int ScreenSeqsCommand::runFastaScreening(map<string, string>& badSeqNames){
	try{
        map<string, string> variables; 
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(fastafile));
        variables["[extension]"] = util.getExtension(fastafile);
		string goodSeqFile = getOutputFileName("fasta", variables);
		outputNames.push_back(goodSeqFile); outputTypes["fasta"].push_back(goodSeqFile);
		
        int numFastaSeqs = createProcesses(goodSeqFile, badAccnosFile, fastafile, badSeqNames);
        
        if (m->getControl_pressed()) { util.mothurRemove(goodSeqFile); return numFastaSeqs; }
		
		return numFastaSeqs;
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "runFastaScreening");
		exit(1);
	}
}
//***************************************************************************************************************/
int ScreenSeqsCommand::screenReports(map<string, string>& badSeqNames){
	try{
        int numFastaSeqs = 0;
        
        //did not provide a summary file, but set a parameter that requires summarizing the fasta file
        //or did provide a summary file, but set maxn parameter so we must summarize the fasta file
        if (((summaryfile == "") && ((util.inUsersGroups("maxambig", optimize)) ||(util.inUsersGroups("maxhomop", optimize)) ||(util.inUsersGroups("maxlength", optimize)) || (util.inUsersGroups("minlength", optimize)) || (util.inUsersGroups("start", optimize)) || (util.inUsersGroups("end", optimize)))) || ((summaryfile != "") && util.inUsersGroups("maxn", optimize))) {  
            getSummary();
        }
        
        if ((summaryfile != "") && ((util.inUsersGroups("maxambig", optimize)) ||(util.inUsersGroups("maxhomop", optimize)) ||(util.inUsersGroups("maxlength", optimize)) || (util.inUsersGroups("minlength", optimize)) || (util.inUsersGroups("start", optimize)) || (util.inUsersGroups("end", optimize)))) { //summarize based on summaryfile
            getSummaryReport();
        }else if ((contigsreport != "") && ((util.inUsersGroups("minoverlap", optimize)) || (util.inUsersGroups("ostart", optimize)) || (util.inUsersGroups("oend", optimize)) || (util.inUsersGroups("mismatches", optimize)))) { //optimize settings based on contigs file
            optimizeContigs();
        }else if ((alignreport != "") && ((util.inUsersGroups("minsim", optimize)) || (util.inUsersGroups("minscore", optimize)) || (util.inUsersGroups("maxinsert", optimize)))) { //optimize settings based on contigs file
            optimizeAlign();
        }
        
        //provided summary file, and did not set maxn so no need to summarize fasta
        if (summaryfile != "")      {   numFastaSeqs = screenSummary(badSeqNames);  }
        //add in any seqs that fail due to contigs report results
        else if (contigsreport != "")    {   numFastaSeqs = screenContigs(badSeqNames);  }
        //add in any seqs that fail due to align report
        else if (alignreport != "")      {   numFastaSeqs = screenAlignReport(badSeqNames);  }
        
        return numFastaSeqs;
    }
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "screenReports");
		exit(1);
	}
}
//***************************************************************************************************************
int ScreenSeqsCommand::screenAlignReport(map<string, string>& badSeqNames){
    try {
        
        map<string, string> variables;
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(alignreport));
        string outSummary =  getOutputFileName("alignreport",variables);
        outputNames.push_back(outSummary); outputTypes["alignreport"].push_back(outSummary);
        
        string name, TemplateName, SearchMethod, AlignmentMethod;
        //QueryName	QueryLength	TemplateName	TemplateLength	SearchMethod	SearchScore	AlignmentMethod	QueryStart	QueryEnd	TemplateStart	TemplateEnd	PairwiseAlignmentLength	GapsInQuery	GapsInTemplate	LongestInsert	SimBtwnQuery&Template
        //checking for minScore, maxInsert, minSim
        int length, TemplateLength,	 QueryStart,	QueryEnd,	TemplateStart,	TemplateEnd,	PairwiseAlignmentLength,	GapsInQuery,	GapsInTemplate,	LongestInsert;
        float SearchScore, SimBtwnQueryTemplate;
        
        ofstream out;
        util.openOutputFile(outSummary, out);
        
        //read summary file
        ifstream in;
        util.openInputFile(alignreport, in);
        out << (util.getline(in)) << endl;   //skip headers
        
        int count = 0;
        
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { in.close(); out.close(); return 0; }
            
            //seqname	start	end	nbases	ambigs	polymer	numSeqs
            in >> name >> length >> TemplateName >> TemplateLength >> SearchMethod >> SearchScore >> AlignmentMethod >> QueryStart >> QueryEnd >> TemplateStart >> TemplateEnd >> PairwiseAlignmentLength >> GapsInQuery >> GapsInTemplate >> LongestInsert >> SimBtwnQueryTemplate; util.gobble(in);
            
            bool goodSeq = 1;		//	innocent until proven guilty
            string trashCode = "";
            if(maxInsert != -1 && maxInsert < LongestInsert)    {	goodSeq = 0; trashCode += "insert|";	}
            if(minScore != -1 && minScore > SearchScore)		{	goodSeq = 0; trashCode += "score|";     }
            if(minSim != -1 && minSim > SimBtwnQueryTemplate)	{	goodSeq = 0; trashCode += "sim|";       }
            
            if(goodSeq == 1){
                out << name << '\t' << length << '\t' << TemplateName  << '\t' << TemplateLength  << '\t' << SearchMethod  << '\t' << SearchScore  << '\t' << AlignmentMethod  << '\t' << QueryStart  << '\t' << QueryEnd  << '\t' << TemplateStart  << '\t' << TemplateEnd  << '\t' << PairwiseAlignmentLength  << '\t' << GapsInQuery  << '\t' << GapsInTemplate  << '\t' << LongestInsert  << '\t' << SimBtwnQueryTemplate << endl;
            }
            else{ badSeqNames[name] = trashCode;  }
            count++;
        }
        in.close();
        out.close();
        
        int oldBadSeqsCount = badSeqNames.size();
        
        int numFastaSeqs = runFastaScreening(badSeqNames);
        
        if (oldBadSeqsCount != badSeqNames.size()) { //more seqs were removed by maxns
            util.renameFile(outSummary, outSummary+".temp");
            
            ofstream out2;
            util.openOutputFile(outSummary, out2);
            
            //read summary file
            ifstream in2;
            util.openInputFile(outSummary+".temp", in2);
            out2 << (util.getline(in2)) << endl;   //skip headers
            
            while (!in2.eof()) {
                
                if (m->getControl_pressed()) { in2.close(); out2.close(); return 0; }
                
                //seqname	start	end	nbases	ambigs	polymer	numSeqs
                in2 >> name >> length >> TemplateName >> TemplateLength >> SearchMethod >> SearchScore >> AlignmentMethod >> QueryStart >> QueryEnd >> TemplateStart >> TemplateEnd >> PairwiseAlignmentLength >> GapsInQuery >> GapsInTemplate >> LongestInsert >> SimBtwnQueryTemplate; util.gobble(in2);
                
                if (badSeqNames.count(name) == 0) { //are you good?
                    out2 << name << '\t' << length << '\t' << TemplateName  << '\t' << TemplateLength  << '\t' << SearchMethod  << '\t' << SearchScore  << '\t' << AlignmentMethod  << '\t' << QueryStart  << '\t' << QueryEnd  << '\t' << TemplateStart  << '\t' << TemplateEnd  << '\t' << PairwiseAlignmentLength  << '\t' << GapsInQuery  << '\t' << GapsInTemplate  << '\t' << LongestInsert  << '\t' << SimBtwnQueryTemplate << endl;
                }
            }
            in2.close();
            out2.close();
            util.mothurRemove(outSummary+".temp");
        }
        
        if (numFastaSeqs != count) {  m->mothurOut("[ERROR]: found " + toString(numFastaSeqs) + " sequences in your fasta file, and " + toString(count) + " sequences in your align report file, quitting.\n"); m->setControl_pressed(true); }
        
        
        return count;
        
        return 0;
        
    }
    catch(exception& e) {
        m->errorOut(e, "ScreenSeqsCommand", "screenAlignReport");
        exit(1);
    }
    
}
//***************************************************************************************************************/
int ScreenSeqsCommand::screenContigs(map<string, string>& badSeqNames){
    try{
        map<string, string> variables;
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(contigsreport));
        string outSummary =  getOutputFileName("contigsreport",variables);
        outputNames.push_back(outSummary); outputTypes["contigsreport"].push_back(outSummary);
        
        string name;
        //Name	Length	Overlap_Length	Overlap_Start	Overlap_End	MisMatches	Num_Ns
        int length, OLength, thisOStart, thisOEnd, numMisMatches, numNs;
        
        ofstream out;
        util.openOutputFile(outSummary, out);
        
        //read summary file
        ifstream in;
        util.openInputFile(contigsreport, in);
        out << (util.getline(in)) << endl;   //skip headers
        
        int count = 0;
        
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { in.close(); out.close(); return 0; }
            
            //seqname	start	end	nbases	ambigs	polymer	numSeqs
            in >> name >> length >> OLength >> thisOStart >> thisOEnd >> numMisMatches >> numNs; util.gobble(in);
            
            bool goodSeq = 1;		//	innocent until proven guilty
            string trashCode = "";
            if(oStart != -1 && oStart < thisOStart)             {	goodSeq = 0;	trashCode += "ostart|";     }
            if(oEnd != -1 && oEnd > thisOEnd)                   {	goodSeq = 0;	trashCode += "oend|";       }
            if(maxN != -1 && maxN <	numNs)                      {	goodSeq = 0;	trashCode += "n|";          }
            if(minOverlap != -1 && minOverlap > OLength)		{	goodSeq = 0;	trashCode += "olength|";    }
            if(mismatches != -1 && mismatches < numMisMatches)	{	goodSeq = 0;	trashCode += "mismatches|"; }
            
            if(goodSeq == 1){
                out << name << '\t' << length  << '\t' << OLength  << '\t' << thisOStart  << '\t' << thisOEnd  << '\t' << numMisMatches  << '\t' << numNs << endl;
            }
            else{ badSeqNames[name] = trashCode; }
            count++;
        }
        in.close();
        out.close();
        
        int oldBadSeqsCount = badSeqNames.size();
        
        int numFastaSeqs = runFastaScreening(badSeqNames);
        
        if (oldBadSeqsCount != badSeqNames.size()) { //more seqs were removed by maxns
            util.renameFile(outSummary, outSummary+".temp");
            
            ofstream out2;
            util.openOutputFile(outSummary, out2);
            
            //read summary file
            ifstream in2;
            util.openInputFile(outSummary+".temp", in2);
            out2 << (util.getline(in2)) << endl;   //skip headers
            
            while (!in2.eof()) {
                
                if (m->getControl_pressed()) { in2.close(); out2.close(); return 0; }
                
                //seqname	start	end	nbases	ambigs	polymer	numSeqs
                in2 >> name >> length >> OLength >> thisOStart >> thisOEnd >> numMisMatches >> numNs; util.gobble(in2);
                
                if (badSeqNames.count(name) == 0) { //are you good?
                    out2 << name << '\t' << length  << '\t' << OLength  << '\t' << thisOStart  << '\t' << thisOEnd  << '\t' << numMisMatches  << '\t' << numNs << endl;
                }
            }
            in2.close();
            out2.close();
            util.mothurRemove(outSummary+".temp");
        }
        
        if (numFastaSeqs != count) {  m->mothurOut("[ERROR]: found " + toString(numFastaSeqs) + " sequences in your fasta file, and " + toString(count) + " sequences in your contigs report file, quitting.\n"); m->setControl_pressed(true); }
        
        
        return count;
        
    }
    catch(exception& e) {
        m->errorOut(e, "ScreenSeqsCommand", "screenContigs");
        exit(1);
    }
}
//***************************************************************************************************************/
int ScreenSeqsCommand::screenSummary(map<string, string>& badSeqNames){
	try{
        map<string, string> variables; 
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(summaryfile));
        string outSummary =  getOutputFileName("summary",variables);
		outputNames.push_back(outSummary); outputTypes["summary"].push_back(outSummary);
        
        string name;
        int start, end, length, ambigs, polymer, numReps;
        
        ofstream out;
        util.openOutputFile(outSummary, out);
                
        //read summary file
        ifstream in;
        util.openInputFile(summaryfile, in);
        out << (util.getline(in)) << endl;   //skip headers
         
		int count = 0;
        
		while (!in.eof()) {
            
            if (m->getControl_pressed()) { in.close(); out.close(); return 0; }
            
            //seqname	start	end	nbases	ambigs	polymer	numSeqs
            in >> name >> start >> end >> length >> ambigs >> polymer >> numReps; util.gobble(in);
            
            bool goodSeq = 1;		//	innocent until proven guilty
            string trashCode = "";
            if(startPos != -1 && startPos < start)			{	goodSeq = 0;	trashCode += "start|"; }
            if(endPos != -1 && endPos > end)				{	goodSeq = 0;	trashCode += "end|"; }
            if(maxAmbig != -1 && maxAmbig <	ambigs)         {	goodSeq = 0;	trashCode += "ambig|"; }
            if(maxHomoP != -1 && maxHomoP < polymer)        {	goodSeq = 0;	trashCode += "homop|"; }
            if(minLength > length)                          {	goodSeq = 0;	trashCode += "<length|"; }
            if(maxLength != -1 && maxLength < length)		{	goodSeq = 0;	trashCode += ">length|"; }
            
            if(goodSeq == 1){
                out << name << '\t' << start  << '\t' << end  << '\t' << length  << '\t' << ambigs  << '\t' << polymer  << '\t' << numReps << endl;	
            }
            else{ badSeqNames[name] = trashCode; }
            count++;
        }
        in.close();
        out.close();
        
        int oldBadSeqsCount = badSeqNames.size();
        
        int numFastaSeqs = runFastaScreening(badSeqNames);
        
        if (oldBadSeqsCount != badSeqNames.size()) { //more seqs were removed by maxns
            util.renameFile(outSummary, outSummary+".temp");
            
            ofstream out2;
            util.openOutputFile(outSummary, out2);
            
            //read summary file
            ifstream in2;
            util.openInputFile(outSummary+".temp", in2);
            out2 << (util.getline(in2)) << endl;   //skip headers
            
            while (!in2.eof()) {
                
                if (m->getControl_pressed()) { in2.close(); out2.close(); return 0; }
                
                //seqname	start	end	nbases	ambigs	polymer	numSeqs
                in2 >> name >> start >> end >> length >> ambigs >> polymer >> numReps; util.gobble(in2);
                
                if (badSeqNames.count(name) == 0) { //are you good?
                    out2 << name << '\t' << start  << '\t' << end  << '\t' << length  << '\t' << ambigs  << '\t' << polymer  << '\t' << numReps << endl;	
                }
            }
            in2.close();
            out2.close();
            util.mothurRemove(outSummary+".temp");
        }
        
        if (numFastaSeqs != count) {  m->mothurOut("[ERROR]: found " + toString(numFastaSeqs) + " sequences in your fasta file, and " + toString(count) + " sequences in your summary file, quitting.\n"); m->setControl_pressed(true); }
        
        
        
        return count;
    }
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "screenSummary");
		exit(1);
	}
}
//***************************************************************************************************************/
int ScreenSeqsCommand::screenFasta(map<string, string>& badSeqNames){
	try{
        if (optimize.size() != 0) {   getSummary();  }
    
        if (m->getControl_pressed()) { return 0; }
        
        int numFastaSeqs = runFastaScreening(badSeqNames);
        
        return numFastaSeqs;
    }
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "screenFasta");
		exit(1);
	}
}	
//***************************************************************************************************************
int ScreenSeqsCommand::getSummaryReport(){
	try {
        Summary sum(processors);
        sum.summarizeFastaSummary(summaryfile);
        
        double criteriaPercentile = criteria;
        double mincriteriaPercentile = (100 - criteria);
        
        for (int i = 0; i < optimize.size(); i++) {
            if (optimize[i] == "start") { startPos = sum.getStart(criteriaPercentile); m->mothurOut("Optimizing start to " + toString(startPos) + "."); m->mothurOutEndLine(); }
            else if (optimize[i] == "end") {   endPos = sum.getEnd(mincriteriaPercentile); m->mothurOut("Optimizing end to " + toString(endPos) + "."); m->mothurOutEndLine();}
            else if (optimize[i] == "maxambig") { maxAmbig = sum.getAmbig(criteriaPercentile); m->mothurOut("Optimizing maxambig to " + toString(maxAmbig) + "."); m->mothurOutEndLine(); }
            else if (optimize[i] == "maxhomop") { maxHomoP = sum.getAmbig(criteriaPercentile); m->mothurOut("Optimizing maxhomop to " + toString(maxHomoP) + "."); m->mothurOutEndLine(); }
            else if (optimize[i] == "minlength") {  minLength = sum.getLength(mincriteriaPercentile); m->mothurOut("Optimizing minlength to " + toString(minLength) + "."); m->mothurOutEndLine(); if (minLength < 0) { m->setControl_pressed(true); } }
            else if (optimize[i] == "maxlength") { maxLength = sum.getLength(criteriaPercentile); m->mothurOut("Optimizing maxlength to " + toString(maxLength) + "."); m->mothurOutEndLine(); }
        }
        
        return 0;
        
    }
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "getSummaryReport");
		exit(1);
	}
}
//***************************************************************************************************************
int ScreenSeqsCommand::optimizeContigs(){
	try{
        Summary sum(processors);
        sum.summarizeContigsSummary(contigsreport);
        
        double criteriaPercentile = criteria;
        double mincriteriaPercentile = (100 - criteria);

        for (int i = 0; i < optimize.size(); i++) {
            if (optimize[i] == "ostart") { oStart = sum.getOStart(criteriaPercentile); m->mothurOut("Optimizing ostart to " + toString(oStart) + "."); m->mothurOutEndLine(); }
            else if (optimize[i] == "oend") {  endPos = sum.getOEnd(mincriteriaPercentile); m->mothurOut("Optimizing oend to " + toString(oEnd) + "."); m->mothurOutEndLine();}
            else if (optimize[i] == "mismatches") { mismatches = sum.getMisMatches(criteriaPercentile); m->mothurOut("Optimizing mismatches to " + toString(mismatches) + "."); m->mothurOutEndLine(); }
            else if (optimize[i] == "maxn") { maxN = sum.getNumNs(criteriaPercentile); m->mothurOut("Optimizing maxn to " + toString(maxN) + "."); m->mothurOutEndLine(); }
            else if (optimize[i] == "minoverlap") {  minOverlap = sum.getOLength(mincriteriaPercentile); m->mothurOut("Optimizing minoverlap to " + toString(minOverlap) + "."); m->mothurOutEndLine(); }
        }
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "optimizeContigs");
		exit(1);
	}
}
//***************************************************************************************************************
int ScreenSeqsCommand::optimizeAlign(){
	try {
        
        Summary sum(processors);
        sum.summarizeAlignSummary(alignreport);
       
        double mincriteriaPercentile = (100 - criteria);
        
        for (int i = 0; i < optimize.size(); i++) {
            if (optimize[i] == "minsim") {  minSim = sum.getSims(mincriteriaPercentile);  m->mothurOut("Optimizing minsim to " + toString(minSim) + "."); m->mothurOutEndLine();}
            else if (optimize[i] == "minscore") {  minScore = sum.getScores(mincriteriaPercentile);  m->mothurOut("Optimizing minscore to " + toString(minScore) + "."); m->mothurOutEndLine(); }
            else if (optimize[i] == "maxinsert") { maxInsert = sum.getNumInserts(mincriteriaPercentile); m->mothurOut("Optimizing maxinsert to " + toString(maxInsert) + "."); m->mothurOutEndLine(); }
        }
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "optimizeAlign");
		exit(1);
	}
}
//***************************************************************************************************************
int ScreenSeqsCommand::getSummary(){
	try {
        Summary sum(processors);
        sum.summarizeFasta(fastafile, "");
        
		//numSeqs is the number of unique seqs, startPosition.size() is the total number of seqs, we want to optimize using all seqs
		double criteriaPercentile = criteria;
        double mincriteriaPercentile = (100 - criteria);
        
        for (int i = 0; i < optimize.size(); i++) {
            if (optimize[i] == "start") { startPos = sum.getStart(criteriaPercentile); m->mothurOut("Optimizing start to " + toString(startPos) + "."); m->mothurOutEndLine(); }
            else if (optimize[i] == "end") {   endPos = sum.getEnd(mincriteriaPercentile); m->mothurOut("Optimizing end to " + toString(endPos) + "."); m->mothurOutEndLine();}
            else if (optimize[i] == "maxambig") { maxAmbig = sum.getAmbig(criteriaPercentile); m->mothurOut("Optimizing maxambig to " + toString(maxAmbig) + "."); m->mothurOutEndLine(); }
            else if (optimize[i] == "maxhomop") { maxHomoP = sum.getAmbig(criteriaPercentile); m->mothurOut("Optimizing maxhomop to " + toString(maxHomoP) + "."); m->mothurOutEndLine(); }
            else if (optimize[i] == "minlength") { minLength = sum.getLength(mincriteriaPercentile); m->mothurOut("Optimizing minlength to " + toString(minLength) + "."); m->mothurOutEndLine(); if (minLength < 0) { m->setControl_pressed(true); } }
            else if (optimize[i] == "maxlength") { maxLength = sum.getLength(criteriaPercentile); m->mothurOut("Optimizing maxlength to " + toString(maxLength) + "."); m->mothurOutEndLine(); }
            else if (optimize[i] == "maxn") { maxN = sum.getNumNs(criteriaPercentile); m->mothurOut("Optimizing maxn to " + toString(maxN) + "."); m->mothurOutEndLine(); }
        }

		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "getSummary");
		exit(1);
	}
}
//**********************************************************************************************************************

void driverScreen(sumScreenData* params){
	try {
		ifstream inFASTA;
		params->util.openInputFile(params->filename, inFASTA);

        inFASTA.seekg(params->start);

        //print header if you are process 0
        if (params->start == 0) { params->util.zapGremlins(inFASTA); params->util.gobble(inFASTA); }

		bool done = false;
		params->count = 0;
        
		while (!done) {
		
			if (params->m->getControl_pressed()) {  break; }
			
			Sequence currSeq(inFASTA); params->util.gobble(inFASTA);
			if (currSeq.getName() != "") {
				bool goodSeq = 1;		//	innocent until proven guilty
                string trashCode = "";
                //have the report files found you bad
                map<string, string>::iterator it = params->badSeqNames.find(currSeq.getName());
                if (it != params->badSeqNames.end()) { goodSeq = 0;  trashCode = it->second; }
                
                if (params->summaryfile == "") { //summaryfile includes these so no need to check again
                    if(params->startPos != -1 && params->startPos < currSeq.getStartPos())			{	goodSeq = 0;	trashCode += "start|";  }
                    if(params->endPos != -1 && params->endPos > currSeq.getEndPos())				{	goodSeq = 0;	trashCode += "end|";    }
                    if(params->maxAmbig != -1 && params->maxAmbig <	currSeq.getAmbigBases())		{	goodSeq = 0;	trashCode += "ambig|";  }
                    if(params->maxHomoP != -1 && params->maxHomoP < currSeq.getLongHomoPolymer())	{	goodSeq = 0;	trashCode += "homop|";  }
                    if(params->minLength > currSeq.getNumBases())                                   {	goodSeq = 0;	trashCode += "<length|";}
                    if(params->maxLength != -1 && params->maxLength < currSeq.getNumBases())		{	goodSeq = 0;	trashCode += ">length|";}
                    
                    if (params->m->getDebug()) { params->m->mothurOut("[DEBUG]: " + currSeq.getName() + "\t" + toString(currSeq.getStartPos()) + "\t" + toString(currSeq.getEndPos()) + "\t" + toString(currSeq.getNumBases()) + "\n"); }
                }
                
                if (params->contigsreport == "") { //contigs report includes this so no need to check again
                    if(params->maxN != -1 && params->maxN < currSeq.getNumNs())                     {	goodSeq = 0;	trashCode += "n|"; }
                }
				
				if(goodSeq == 1){
					currSeq.printSequence(params->outputWriter);
				}else{
					string badAccnos = currSeq.getName() + '\t' + trashCode.substr(0, trashCode.length()-1) + '\n';
                    params->accnosWriter->write(badAccnos);
					params->badSeqNames[currSeq.getName()] = trashCode;
				}
                params->count++;
			}
			
			#if defined NON_WINDOWS
				unsigned long long pos = inFASTA.tellg();
				if ((pos == -1) || (pos >= params->end)) { break; }
			#else
				if (params->end == params->count) { break; }
			#endif
			
			//report progress
			if((params->count) % 1000 == 0){	params->m->mothurOutJustToScreen(toString(params->count)+"\n"); 		}
		}
		//report progress
		if((params->count) % 1000 != 0){	params->m->mothurOutJustToScreen(toString(params->count)+"\n"); 	}
		
		inFASTA.close();
	}
	catch(exception& e) {
		params->m->errorOut(e, "ScreenSeqsCommand", "driverScreen");
		exit(1);
	}
}
/**************************************************************************************************/

int ScreenSeqsCommand::createProcesses(string goodFileName, string badAccnos, string filename, map<string, string>& badSeqNames) {
	try {
        
        vector<linePair> lines;
        vector<unsigned long long> positions;
#if defined NON_WINDOWS
        positions = util.divideFile(fastafile, processors);
        for (int i = 0; i < (positions.size()-1); i++) { lines.push_back(linePair(positions[i], positions[(i+1)])); }
#else
        
        long long numFastaSeqs = 0;
        positions = util.setFilePosFasta(fastafile, numFastaSeqs);
        if (numFastaSeqs < processors) { processors = numFastaSeqs; }
        
        //figure out how many sequences you have to process
        int numSeqsPerProcessor = numFastaSeqs / processors;
        for (int i = 0; i < processors; i++) {
            int startIndex =  i * numSeqsPerProcessor;
            if(i == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor; 	}
            lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
        }
        
#endif
        
        //create array of worker threads
        vector<thread*> workerThreads;
        vector<sumScreenData*> data;
        
        long long num = 0;
        
        time_t start, end;
        time(&start);
        
        auto synchronizedOutputFile = std::make_shared<SynchronizedOutputFile>(goodFileName);
        auto synchronizedAccnosFile = std::make_shared<SynchronizedOutputFile>(badAccnos);
        
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            
            OutputWriter* outputThreadWriter = new OutputWriter(synchronizedOutputFile);
            OutputWriter* accnosThreadWriter = new OutputWriter(synchronizedAccnosFile);
            
            sumScreenData* dataBundle = new sumScreenData(startPos, endPos, maxAmbig, maxHomoP, minLength, maxLength, maxN, badSeqNames, filename, summaryfile, contigsreport, lines[i+1].start, lines[i+1].end,outputThreadWriter, accnosThreadWriter);
            
            data.push_back(dataBundle);
            workerThreads.push_back(new thread(driverScreen, dataBundle));
        }
        
        OutputWriter* outputThreadWriter = new OutputWriter(synchronizedOutputFile);
        OutputWriter* accnosThreadWriter = new OutputWriter(synchronizedAccnosFile);

        sumScreenData* dataBundle = new sumScreenData(startPos, endPos, maxAmbig, maxHomoP, minLength, maxLength, maxN, badSeqNames, filename, summaryfile, contigsreport, lines[0].start, lines[0].end,outputThreadWriter, accnosThreadWriter);
        driverScreen(dataBundle);
        num = dataBundle->count;
        for (map<string, string>::iterator it = dataBundle->badSeqNames.begin(); it != dataBundle->badSeqNames.end(); it++) {	badSeqNames[it->first] = it->second;       }
  
        
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            num += data[i]->count;
            
            for (map<string, string>::iterator it = data[i]->badSeqNames.begin(); it != data[i]->badSeqNames.end(); it++) {	badSeqNames[it->first] = it->second;       }
            
            delete data[i]->outputWriter;
            delete data[i]->accnosWriter;
            delete data[i];
            delete workerThreads[i];
        }
        long long numRemoved = badSeqNames.size();
        
        time(&end);
        m->mothurOut("\nIt took " + toString(difftime(end, start)) + " secs to screen " + toString(num) + " sequences, removed " + toString(numRemoved) + ".\n\n");
        delete outputThreadWriter; delete accnosThreadWriter;
        delete dataBundle;
        return num;
        
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "createProcesses");
		exit(1);
	}
}

//***************************************************************************************************************


