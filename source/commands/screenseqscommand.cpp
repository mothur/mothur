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
		CommandParameter pminlength("minlength", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pminlength);
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
		helpString += "The minlength parameter allows you to set and minimum sequence length. \n";
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
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFasta).\n";
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
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
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
			string inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("fasta");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fasta"] = inputDir + it->second;		}
				}
				
				it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
				
				it = parameters.find("name");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["name"] = inputDir + it->second;		}
				}
				
				it = parameters.find("alignreport");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["alignreport"] = inputDir + it->second;		}
				}
                
                it = parameters.find("contigsreport");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["contigsreport"] = inputDir + it->second;		}
				}
                
                it = parameters.find("summary");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["summary"] = inputDir + it->second;		}
				}
				
				it = parameters.find("qfile");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["qfile"] = inputDir + it->second;		}
				}
				
				it = parameters.find("taxonomy");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["taxonomy"] = inputDir + it->second;		}
				}
                
                it = parameters.find("count");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["count"] = inputDir + it->second;		}
				}
			}

			//check for required parameters
			fastafile = validParameter.validFile(parameters, "fasta", true);
			if (fastafile == "not found") { 			
				fastafile = m->getFastaFile(); 
				if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}
			else if (fastafile == "not open") { abort = true; }
			else { m->setFastaFile(fastafile); }
	
			groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") { abort = true; }	
			else if (groupfile == "not found") { groupfile = ""; }
			else { m->setGroupFile(groupfile); }
			
			qualfile = validParameter.validFile(parameters, "qfile", true);
			if (qualfile == "not open") { abort = true; }	
			else if (qualfile == "not found") { qualfile = ""; }
			else { m->setQualFile(qualfile); }
			
			namefile = validParameter.validFile(parameters, "name", true);
			if (namefile == "not open") { namefile = ""; abort = true; }
			else if (namefile == "not found") { namefile = ""; }	
			else { m->setNameFile(namefile); }
			
            countfile = validParameter.validFile(parameters, "count", true);
			if (countfile == "not open") { countfile = ""; abort = true; }
			else if (countfile == "not found") { countfile = "";  }	
			else { m->setCountTableFile(countfile); }
            
            contigsreport = validParameter.validFile(parameters, "contigsreport", true);
			if (contigsreport == "not open") { contigsreport = ""; abort = true; }
			else if (contigsreport == "not found") { contigsreport = "";  }	
            
            summaryfile = validParameter.validFile(parameters, "summary", true);
			if (summaryfile == "not open") { summaryfile = ""; abort = true; }
			else if (summaryfile == "not found") { summaryfile = "";  }
            else { m->setSummaryFile(summaryfile); }
            
            if ((namefile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: name or count."); m->mothurOutEndLine(); abort = true;
            }
			
            if ((groupfile != "") && (countfile != "")) {
                m->mothurOut("[ERROR]: you may only use one of the following: group or count."); m->mothurOutEndLine(); abort=true;
            }
            
			alignreport = validParameter.validFile(parameters, "alignreport", true);
			if (alignreport == "not open") { abort = true; }
			else if (alignreport == "not found") { alignreport = ""; }
			
			taxonomy = validParameter.validFile(parameters, "taxonomy", true);
			if (taxonomy == "not open") { abort = true; }
			else if (taxonomy == "not found") { taxonomy = ""; }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += m->hasPath(fastafile); //if user entered a file with a path then preserve it	
			}

			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			string temp;
			temp = validParameter.validFile(parameters, "start", false);		if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, startPos); 
		
			temp = validParameter.validFile(parameters, "end", false);			if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, endPos);  

			temp = validParameter.validFile(parameters, "maxambig", false);		if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, maxAmbig);  

			temp = validParameter.validFile(parameters, "maxhomop", false);		if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, maxHomoP);  

			temp = validParameter.validFile(parameters, "minlength", false);	if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, minLength); 
			
			temp = validParameter.validFile(parameters, "maxlength", false);	if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, maxLength); 
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors);
			
            temp = validParameter.validFile(parameters, "minoverlap", false);	if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, minOverlap); 
            
            temp = validParameter.validFile(parameters, "ostart", false);	if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, oStart); 
            
            temp = validParameter.validFile(parameters, "oend", false);	if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, oEnd); 
            
            temp = validParameter.validFile(parameters, "mismatches", false);	if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, mismatches); 
            
            temp = validParameter.validFile(parameters, "maxn", false);	if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, maxN); 
            
            temp = validParameter.validFile(parameters, "minscore", false);	if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, minScore); 
            
            temp = validParameter.validFile(parameters, "maxinsert", false);	if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, maxInsert); 
            
            temp = validParameter.validFile(parameters, "minsim", false);	if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, minSim); 
            
			temp = validParameter.validFile(parameters, "optimize", false);	//optimizing trumps the optimized values original value
			if (temp == "not found"){	temp = "none";		}
			m->splitAtDash(temp, optimize);		
            
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
			
			temp = validParameter.validFile(parameters, "criteria", false);	if (temp == "not found"){	temp = "90";				}
			m->mothurConvert(temp, criteria); 
			
			if (countfile == "") { 
                if (namefile == "") {
                    vector<string> files; files.push_back(fastafile);
                    parser.getNameFile(files);
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
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
        map<string, string> badSeqNames;
        int start = time(NULL);
        int numFastaSeqs = 0;
        
        if ((contigsreport == "") && (summaryfile == "") && (alignreport == "")) {   numFastaSeqs = screenFasta(badSeqNames);  }
        else {   numFastaSeqs = screenReports(badSeqNames);   }
		
        if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]); } return 0; }
        
        #ifdef USE_MPI
            int pid;
            MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
        
            if (pid == 0) { //only one process should fix files
        #endif	
                
		if(namefile != "" && groupfile != "")	{	
			screenNameGroupFile(badSeqNames);	
			if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]);  } return 0; }
		}else if(namefile != "")	{	
			screenNameGroupFile(badSeqNames);
			if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]);  } return 0; }	
		}else if(groupfile != "")				{	screenGroupFile(badSeqNames);		}	// this screens just the group
		else if (countfile != "") {     screenCountFile(badSeqNames);		}
            
                
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]);  } return 0; }

		if(qualfile != "")						{	screenQual(badSeqNames);			}
		if(taxonomy != "")						{	screenTaxonomy(badSeqNames);		}
		
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) { m->mothurRemove(outputNames[i]);  } return 0; }
		
		#ifdef USE_MPI
			}
		#endif

		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) { m->mothurOut(outputNames[i]); m->mothurOutEndLine(); }
		m->mothurOutEndLine();
		m->mothurOutEndLine();
		
		//set fasta file as new current fastafile
		string current = "";
		itTypes = outputTypes.find("fasta");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setFastaFile(current); }
		}
		
		itTypes = outputTypes.find("name");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setNameFile(current); }
		}
		
		itTypes = outputTypes.find("group");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setGroupFile(current); }
		}
		
		itTypes = outputTypes.find("qfile");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setQualFile(current); }
		}
		
		itTypes = outputTypes.find("taxonomy");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setTaxonomyFile(current); }
		}
        
        itTypes = outputTypes.find("count");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { current = (itTypes->second)[0]; m->setCountTableFile(current); }
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
        int numFastaSeqs = 0;
        map<string, string> variables; 
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(fastafile));
        string badAccnosFile =  getOutputFileName("accnos",variables);
        variables["[extension]"] = m->getExtension(fastafile);
		string goodSeqFile = getOutputFileName("fasta", variables);
		outputNames.push_back(goodSeqFile); outputTypes["fasta"].push_back(goodSeqFile);
		outputNames.push_back(badAccnosFile); outputTypes["accnos"].push_back(badAccnosFile);
        
#ifdef USE_MPI	
        int pid, numSeqsPerProcessor; 
        int tag = 2001;
        vector<unsigned long long> MPIPos;
        
        MPI_Status status; 
        MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are
        MPI_Comm_size(MPI_COMM_WORLD, &processors); 
        
        MPI_File inMPI;
        MPI_File outMPIGood;
        MPI_File outMPIBadAccnos;
        
        int outMode=MPI_MODE_CREATE|MPI_MODE_WRONLY; 
        int inMode=MPI_MODE_RDONLY; 
        
        char outGoodFilename[1024];
        strcpy(outGoodFilename, goodSeqFile.c_str());
        
        char outBadAccnosFilename[1024];
        strcpy(outBadAccnosFilename, badAccnosFile.c_str());
        
        char inFileName[1024];
        strcpy(inFileName, fastafile.c_str());
        
        MPI_File_open(MPI_COMM_WORLD, inFileName, inMode, MPI_INFO_NULL, &inMPI);  //comm, filename, mode, info, filepointer
        MPI_File_open(MPI_COMM_WORLD, outGoodFilename, outMode, MPI_INFO_NULL, &outMPIGood);
        MPI_File_open(MPI_COMM_WORLD, outBadAccnosFilename, outMode, MPI_INFO_NULL, &outMPIBadAccnos);
        
        if (m->control_pressed) { MPI_File_close(&inMPI);  MPI_File_close(&outMPIGood); MPI_File_close(&outMPIBadAccnos); return 0; }
        
        if (pid == 0) { //you are the root process 
            
            MPIPos = m->setFilePosFasta(fastafile, numFastaSeqs); //fills MPIPos, returns numSeqs
            
            //send file positions to all processes
            for(int i = 1; i < processors; i++) { 
                MPI_Send(&numFastaSeqs, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
                MPI_Send(&MPIPos[0], (numFastaSeqs+1), MPI_LONG, i, tag, MPI_COMM_WORLD);
            }
            
            //figure out how many sequences you have to align
            numSeqsPerProcessor = numFastaSeqs / processors;
            int startIndex =  pid * numSeqsPerProcessor;
            if(pid == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - pid * numSeqsPerProcessor; 	}
            
            //align your part
            driverMPI(startIndex, numSeqsPerProcessor, inMPI, outMPIGood, outMPIBadAccnos, MPIPos, badSeqNames);
            
            if (m->control_pressed) { MPI_File_close(&inMPI);  MPI_File_close(&outMPIGood);  MPI_File_close(&outMPIBadAccnos);  return 0; }
            
            for (int i = 1; i < processors; i++) {
                //get bad lists
                int badSize;
                MPI_Recv(&badSize, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
            }
        }else{ //you are a child process
            MPI_Recv(&numFastaSeqs, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
            MPIPos.resize(numFastaSeqs+1);
            MPI_Recv(&MPIPos[0], (numFastaSeqs+1), MPI_LONG, 0, tag, MPI_COMM_WORLD, &status);
            
            //figure out how many sequences you have to align
            numSeqsPerProcessor = numFastaSeqs / processors;
            int startIndex =  pid * numSeqsPerProcessor;
            if(pid == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - pid * numSeqsPerProcessor; 	}
            
            //align your part
            driverMPI(startIndex, numSeqsPerProcessor, inMPI, outMPIGood, outMPIBadAccnos, MPIPos, badSeqNames);
            
            if (m->control_pressed) { MPI_File_close(&inMPI);  MPI_File_close(&outMPIGood);  MPI_File_close(&outMPIBadAccnos); return 0; }
            
            //send bad list	
            int badSize = badSeqNames.size();
            MPI_Send(&badSize, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
        }
        
        //close files 
        MPI_File_close(&inMPI);
        MPI_File_close(&outMPIGood);
        MPI_File_close(&outMPIBadAccnos);
        MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
        
#else
        if(processors == 1){ numFastaSeqs = driver(lines[0], goodSeqFile, badAccnosFile, fastafile, badSeqNames);	}	
        else{ numFastaSeqs = createProcesses(goodSeqFile, badAccnosFile, fastafile, badSeqNames); }
        
        if (m->control_pressed) { m->mothurRemove(goodSeqFile); return numFastaSeqs; }
#endif		
        
#ifdef USE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
        
        if (pid == 0) { //only one process should fix files
			
            //read accnos file with all names in it, process 0 just has its names
            MPI_File inMPIAccnos;
            MPI_Offset size;
			
            char inFileName[1024];
            strcpy(inFileName, badAccnosFile.c_str());
			
            MPI_File_open(MPI_COMM_SELF, inFileName, inMode, MPI_INFO_NULL, &inMPIAccnos);  //comm, filename, mode, info, filepointer
            MPI_File_get_size(inMPIAccnos, &size);
			
            char* buffer = new char[size];
            MPI_File_read(inMPIAccnos, buffer, size, MPI_CHAR, &status);
			
            string tempBuf = buffer;
            if (tempBuf.length() > size) { tempBuf = tempBuf.substr(0, size);  }
            istringstream iss (tempBuf,istringstream::in);
            
            delete buffer;
            MPI_File_close(&inMPIAccnos);
            
            badSeqNames.clear();
            string tempName, trashCode;
            while (!iss.eof()) {
                iss >> tempName >> trashCode; m->gobble(iss);
                badSeqNames[tempName] = trashCode;
            }
        }
#endif
        
        
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
        bool summarizedFasta = false;
        
        //did not provide a summary file, but set a parameter that requires summarizing the fasta file
        //or did provide a summary file, but set maxn parameter so we must summarize the fasta file 
        vector<unsigned long long> positions;
        if (((summaryfile == "") && ((m->inUsersGroups("maxambig", optimize)) ||(m->inUsersGroups("maxhomop", optimize)) ||(m->inUsersGroups("maxlength", optimize)) || (m->inUsersGroups("minlength", optimize)) || (m->inUsersGroups("start", optimize)) || (m->inUsersGroups("end", optimize)))) || ((summaryfile != "") && m->inUsersGroups("maxn", optimize))) {  
            //use the namefile to optimize correctly
            if (namefile != "") { nameMap = m->readNames(namefile); }
            else if (countfile != "") {
                CountTable ct;
                ct.readTable(countfile, true, false);
                nameMap = ct.getNameMap();
            }
            getSummary(positions); 
            summarizedFasta = true;
        } else {
            #if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
                positions = m->divideFile(fastafile, processors);
                for (int i = 0; i < (positions.size()-1); i++) { lines.push_back(linePair(positions[i], positions[(i+1)])); }
            #else 
                if(processors == 1){ lines.push_back(linePair(0, 1000));  }
                else {
                    int numFastaSeqs = 0;
                    positions = m->setFilePosFasta(fastafile, numFastaSeqs); 
                    if (positions.size() < processors) { processors = positions.size(); }
                
                    //figure out how many sequences you have to process
                    int numSeqsPerProcessor = numFastaSeqs / processors;
                    for (int i = 0; i < processors; i++) {
                        int startIndex =  i * numSeqsPerProcessor;
                        if(i == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor; 	}
                        lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
                    }
                }
            #endif
        }
        
        if ((summaryfile != "") && ((m->inUsersGroups("maxambig", optimize)) ||(m->inUsersGroups("maxhomop", optimize)) ||(m->inUsersGroups("maxlength", optimize)) || (m->inUsersGroups("minlength", optimize)) || (m->inUsersGroups("start", optimize)) || (m->inUsersGroups("end", optimize))) && !summarizedFasta) { //summarize based on summaryfile
            if (namefile != "") { nameMap = m->readNames(namefile); }
            else if (countfile != "") {
                CountTable ct;
                ct.readTable(countfile, true, false);
                nameMap = ct.getNameMap();
            }
            getSummaryReport();
        }else if ((contigsreport != "") && ((m->inUsersGroups("minoverlap", optimize)) || (m->inUsersGroups("ostart", optimize)) || (m->inUsersGroups("oend", optimize)) || (m->inUsersGroups("mismatches", optimize)))) { //optimize settings based on contigs file
            optimizeContigs();
        }else if ((alignreport != "") && ((m->inUsersGroups("minsim", optimize)) || (m->inUsersGroups("minscore", optimize)) || (m->inUsersGroups("maxinsert", optimize)))) { //optimize settings based on contigs file
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
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(alignreport));
        string outSummary =  getOutputFileName("alignreport",variables);
		outputNames.push_back(outSummary); outputTypes["alignreport"].push_back(outSummary);
        
        string name, TemplateName, SearchMethod, AlignmentMethod;
        //QueryName	QueryLength	TemplateName	TemplateLength	SearchMethod	SearchScore	AlignmentMethod	QueryStart	QueryEnd	TemplateStart	TemplateEnd	PairwiseAlignmentLength	GapsInQuery	GapsInTemplate	LongestInsert	SimBtwnQuery&Template
        //checking for minScore, maxInsert, minSim
        int length, TemplateLength,	 QueryStart,	QueryEnd,	TemplateStart,	TemplateEnd,	PairwiseAlignmentLength,	GapsInQuery,	GapsInTemplate,	LongestInsert;
        float SearchScore, SimBtwnQueryTemplate;
        
        ofstream out;
        m->openOutputFile(outSummary, out);
        
        //read summary file
        ifstream in;
        m->openInputFile(alignreport, in);
        out << (m->getline(in)) << endl;   //skip headers
        
		int count = 0;
        
		while (!in.eof()) {
            
            if (m->control_pressed) { in.close(); out.close(); return 0; }
            
            //seqname	start	end	nbases	ambigs	polymer	numSeqs
            in >> name >> length >> TemplateName >> TemplateLength >> SearchMethod >> SearchScore >> AlignmentMethod >> QueryStart >> QueryEnd >> TemplateStart >> TemplateEnd >> PairwiseAlignmentLength >> GapsInQuery >> GapsInTemplate >> LongestInsert >> SimBtwnQueryTemplate; m->gobble(in);

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
            m->renameFile(outSummary, outSummary+".temp");
            
            ofstream out2;
            m->openOutputFile(outSummary, out2);
            
            //read summary file
            ifstream in2;
            m->openInputFile(outSummary+".temp", in2);
            out2 << (m->getline(in2)) << endl;   //skip headers
            
            while (!in2.eof()) {
                
                if (m->control_pressed) { in2.close(); out2.close(); return 0; }
                
                //seqname	start	end	nbases	ambigs	polymer	numSeqs
                in2 >> name >> length >> TemplateName >> TemplateLength >> SearchMethod >> SearchScore >> AlignmentMethod >> QueryStart >> QueryEnd >> TemplateStart >> TemplateEnd >> PairwiseAlignmentLength >> GapsInQuery >> GapsInTemplate >> LongestInsert >> SimBtwnQueryTemplate; m->gobble(in2);
                
                if (badSeqNames.count(name) == 0) { //are you good?
                    out2 << name << '\t' << length << '\t' << TemplateName  << '\t' << TemplateLength  << '\t' << SearchMethod  << '\t' << SearchScore  << '\t' << AlignmentMethod  << '\t' << QueryStart  << '\t' << QueryEnd  << '\t' << TemplateStart  << '\t' << TemplateEnd  << '\t' << PairwiseAlignmentLength  << '\t' << GapsInQuery  << '\t' << GapsInTemplate  << '\t' << LongestInsert  << '\t' << SimBtwnQueryTemplate << endl;		
                }
            }
            in2.close();
            out2.close();
            m->mothurRemove(outSummary+".temp");
        }
        
        if (numFastaSeqs != count) {  m->mothurOut("[ERROR]: found " + toString(numFastaSeqs) + " sequences in your fasta file, and " + toString(count) + " sequences in your align report file, quitting.\n"); m->control_pressed = true; }
        
        
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
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(contigsreport));
        string outSummary =  getOutputFileName("contigsreport",variables);
		outputNames.push_back(outSummary); outputTypes["contigsreport"].push_back(outSummary);
        
        string name;
        //Name	Length	Overlap_Length	Overlap_Start	Overlap_End	MisMatches	Num_Ns
        int length, OLength, thisOStart, thisOEnd, numMisMatches, numNs;
        
        ofstream out;
        m->openOutputFile(outSummary, out);
        
        //read summary file
        ifstream in;
        m->openInputFile(contigsreport, in);
        out << (m->getline(in)) << endl;   //skip headers
        
		int count = 0;
        
		while (!in.eof()) {
            
            if (m->control_pressed) { in.close(); out.close(); return 0; }
            
            //seqname	start	end	nbases	ambigs	polymer	numSeqs
            in >> name >> length >> OLength >> thisOStart >> thisOEnd >> numMisMatches >> numNs; m->gobble(in);
            
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
            m->renameFile(outSummary, outSummary+".temp");
            
            ofstream out2;
            m->openOutputFile(outSummary, out2);
            
            //read summary file
            ifstream in2;
            m->openInputFile(outSummary+".temp", in2);
            out2 << (m->getline(in2)) << endl;   //skip headers
            
            while (!in2.eof()) {
                
                if (m->control_pressed) { in2.close(); out2.close(); return 0; }
                
                //seqname	start	end	nbases	ambigs	polymer	numSeqs
                in2 >> name >> length >> OLength >> thisOStart >> thisOEnd >> numMisMatches >> numNs; m->gobble(in2);
                
                if (badSeqNames.count(name) == 0) { //are you good?
                    out2 << name << '\t' << length  << '\t' << OLength  << '\t' << thisOStart  << '\t' << thisOEnd  << '\t' << numMisMatches  << '\t' << numNs << endl;		
                }
            }
            in2.close();
            out2.close();
            m->mothurRemove(outSummary+".temp");
        }
        
        if (numFastaSeqs != count) {  m->mothurOut("[ERROR]: found " + toString(numFastaSeqs) + " sequences in your fasta file, and " + toString(count) + " sequences in your contigs report file, quitting.\n"); m->control_pressed = true; }
        
        
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
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(summaryfile));
        string outSummary =  getOutputFileName("summary",variables);
		outputNames.push_back(outSummary); outputTypes["summary"].push_back(outSummary);
        
        string name;
        int start, end, length, ambigs, polymer, numReps;
        
        ofstream out;
        m->openOutputFile(outSummary, out);
                
        //read summary file
        ifstream in;
        m->openInputFile(summaryfile, in);
        out << (m->getline(in)) << endl;   //skip headers
         
		int count = 0;
        
		while (!in.eof()) {
            
            if (m->control_pressed) { in.close(); out.close(); return 0; }
            
            //seqname	start	end	nbases	ambigs	polymer	numSeqs
            in >> name >> start >> end >> length >> ambigs >> polymer >> numReps; m->gobble(in);
            
            bool goodSeq = 1;		//	innocent until proven guilty
            string trashCode = "";
            if(startPos != -1 && startPos < start)			{	goodSeq = 0;	trashCode += "start|"; }
            if(endPos != -1 && endPos > end)				{	goodSeq = 0;	trashCode += "end|"; }
            if(maxAmbig != -1 && maxAmbig <	ambigs)         {	goodSeq = 0;	trashCode += "ambig|"; }
            if(maxHomoP != -1 && maxHomoP < polymer)        {	goodSeq = 0;	trashCode += "homop|"; }
            if(minLength != -1 && minLength > length)		{	goodSeq = 0;	trashCode += "<length|"; }
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
            m->renameFile(outSummary, outSummary+".temp");
            
            ofstream out2;
            m->openOutputFile(outSummary, out2);
            
            //read summary file
            ifstream in2;
            m->openInputFile(outSummary+".temp", in2);
            out2 << (m->getline(in2)) << endl;   //skip headers
            
            while (!in2.eof()) {
                
                if (m->control_pressed) { in2.close(); out2.close(); return 0; }
                
                //seqname	start	end	nbases	ambigs	polymer	numSeqs
                in2 >> name >> start >> end >> length >> ambigs >> polymer >> numReps; m->gobble(in2);
                
                if (badSeqNames.count(name) == 0) { //are you good?
                    out2 << name << '\t' << start  << '\t' << end  << '\t' << length  << '\t' << ambigs  << '\t' << polymer  << '\t' << numReps << endl;	
                }
            }
            in2.close();
            out2.close();
            m->mothurRemove(outSummary+".temp");
        }
        
        if (numFastaSeqs != count) {  m->mothurOut("[ERROR]: found " + toString(numFastaSeqs) + " sequences in your fasta file, and " + toString(count) + " sequences in your summary file, quitting.\n"); m->control_pressed = true; }
        
        
        
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
        
        
        //if the user want to optimize we need to know the 90% mark
		vector<unsigned long long> positions;
		if (optimize.size() != 0) {  //get summary is paralellized so we need to divideFile, no need to do this step twice so I moved it here
			//use the namefile to optimize correctly
			if (namefile != "") { nameMap = m->readNames(namefile); }
            else if (countfile != "") {
                CountTable ct;
                ct.readTable(countfile, true, false);
                nameMap = ct.getNameMap();
            }
			getSummary(positions); 
		}else { 
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
            positions = m->divideFile(fastafile, processors);
            for (int i = 0; i < (positions.size()-1); i++) { lines.push_back(linePair(positions[i], positions[(i+1)])); }
#else 
            if(processors == 1){ lines.push_back(linePair(0, 1000));  }
            else {
                int numFastaSeqs = 0;
                positions = m->setFilePosFasta(fastafile, numFastaSeqs); 
                if (positions.size() < processors) { processors = positions.size(); }
                
                //figure out how many sequences you have to process
                int numSeqsPerProcessor = numFastaSeqs / processors;
                for (int i = 0; i < processors; i++) {
                    int startIndex =  i * numSeqsPerProcessor;
                    if(i == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor; 	}
                    lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
                }
            }
#endif
		}
        
        if (m->control_pressed) { return 0; }
        
        int numFastaSeqs = runFastaScreening(badSeqNames);
        
        return numFastaSeqs;
        
    }
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "screenFasta");
		exit(1);
	}
}	
//***************************************************************************************************************

int ScreenSeqsCommand::screenNameGroupFile(map<string, string> badSeqNames){
	try {
		ifstream inputNames;
		m->openInputFile(namefile, inputNames);
		map<string, string> badSeqGroups;
		string seqName, seqList, group;
		map<string, string>::iterator it;
        map<string, string> variables; 
		variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(namefile));
        variables["[extension]"] = m->getExtension(namefile);
		string goodNameFile = getOutputFileName("name", variables);
		outputNames.push_back(goodNameFile);  outputTypes["name"].push_back(goodNameFile);
		
		ofstream goodNameOut;	m->openOutputFile(goodNameFile, goodNameOut);
		
		while(!inputNames.eof()){
			if (m->control_pressed) { goodNameOut.close();  inputNames.close(); m->mothurRemove(goodNameFile);  return 0; }

			inputNames >> seqName; m->gobble(inputNames); inputNames >> seqList;
			it = badSeqNames.find(seqName);
				
			if(it != badSeqNames.end()){
				if(namefile != ""){
					int start = 0;
					for(int i=0;i<seqList.length();i++){
						if(seqList[i] == ','){
							badSeqGroups[seqList.substr(start,i-start)] = it->second;
							start = i+1;
						}					
					}
					badSeqGroups[seqList.substr(start,seqList.length()-start)] = it->second;
				}
                badSeqNames.erase(it);
			}
			else{
				goodNameOut << seqName << '\t' << seqList << endl;
			}
			m->gobble(inputNames);
		}
		inputNames.close();
		goodNameOut.close();
	
		//we were unable to remove some of the bad sequences
		if (badSeqNames.size() != 0) {
			for (it = badSeqNames.begin(); it != badSeqNames.end(); it++) {  
				m->mothurOut("Your namefile does not include the sequence " + it->first + " please correct."); 
				m->mothurOutEndLine();
			}
		}

		if(groupfile != ""){
			
			ifstream inputGroups;
			m->openInputFile(groupfile, inputGroups);
            variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(groupfile));
            variables["[extension]"] = m->getExtension(groupfile);
            string goodGroupFile = getOutputFileName("group", variables);
			
			outputNames.push_back(goodGroupFile);   outputTypes["group"].push_back(goodGroupFile);
			
			ofstream goodGroupOut;	m->openOutputFile(goodGroupFile, goodGroupOut);
			
			while(!inputGroups.eof()){
				if (m->control_pressed) { goodGroupOut.close(); inputGroups.close(); m->mothurRemove(goodNameFile);  m->mothurRemove(goodGroupFile); return 0; }

				inputGroups >> seqName; m->gobble(inputGroups); inputGroups >> group;
				
				it = badSeqGroups.find(seqName);
				
				if(it != badSeqGroups.end()){
					badSeqGroups.erase(it);
				}
				else{
					goodGroupOut << seqName << '\t' << group << endl;
				}
				m->gobble(inputGroups);
			}
			inputGroups.close();
			goodGroupOut.close();
			
			//we were unable to remove some of the bad sequences
			if (badSeqGroups.size() != 0) {
				for (it = badSeqGroups.begin(); it != badSeqGroups.end(); it++) {  
					m->mothurOut("Your groupfile does not include the sequence " + it->first + " please correct."); 
					m->mothurOutEndLine();
				}
			}
		}
		
		
		return 0;
	
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "screenNameGroupFile");
		exit(1);
	}
}
//***************************************************************************************************************
int ScreenSeqsCommand::getSummaryReport(){
	try {
		
		vector<int> startPosition;
		vector<int> endPosition;
		vector<int> seqLength;
		vector<int> ambigBases;
		vector<int> longHomoPolymer;
        
#ifdef USE_MPI
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
		
		if (pid == 0) { 
#endif
            
            
            //read summary file
            ifstream in;
            m->openInputFile(summaryfile, in);
            m->getline(in);
            
            string name;
            int start, end, length, ambigs, polymer, numReps;
            
            while (!in.eof()) {
                
                if (m->control_pressed) { in.close(); return 0; }
                
                //seqname	start	end	nbases	ambigs	polymer	numSeqs
                in >> name >> start >> end >> length >> ambigs >> polymer >> numReps; m->gobble(in);
                
                int num = 1;
				if ((namefile != "") || (countfile !="")) {
					//make sure this sequence is in the namefile, else error 
					map<string, int>::iterator it = nameMap.find(name);
					
					if (it == nameMap.end()) { m->mothurOut("[ERROR]: " + name + " is not in your namefile, please correct."); m->mothurOutEndLine(); m->control_pressed = true; }
					else { num = it->second; }
				}
				
				//for each sequence this sequence represents
				for (int i = 0; i < num; i++) {
					startPosition.push_back(start);
					endPosition.push_back(end);
					seqLength.push_back(length);
					ambigBases.push_back(ambigs);
					longHomoPolymer.push_back(polymer);
				}
               
            }
            in.close();

        sort(startPosition.begin(), startPosition.end());
		sort(endPosition.begin(), endPosition.end());
		sort(seqLength.begin(), seqLength.end());
		sort(ambigBases.begin(), ambigBases.end());
		sort(longHomoPolymer.begin(), longHomoPolymer.end());
		
		//numSeqs is the number of unique seqs, startPosition.size() is the total number of seqs, we want to optimize using all seqs
		int criteriaPercentile	= int(startPosition.size() * (criteria / (float) 100));
		
		for (int i = 0; i < optimize.size(); i++) {
			if (optimize[i] == "start") { startPos = startPosition[criteriaPercentile]; m->mothurOut("Optimizing start to " + toString(startPos) + "."); m->mothurOutEndLine(); }
			else if (optimize[i] == "end") { int endcriteriaPercentile = int(endPosition.size() * ((100 - criteria) / (float) 100));  endPos = endPosition[endcriteriaPercentile]; m->mothurOut("Optimizing end to " + toString(endPos) + "."); m->mothurOutEndLine();}
			else if (optimize[i] == "maxambig") { maxAmbig = ambigBases[criteriaPercentile]; m->mothurOut("Optimizing maxambig to " + toString(maxAmbig) + "."); m->mothurOutEndLine(); }
			else if (optimize[i] == "maxhomop") { maxHomoP = longHomoPolymer[criteriaPercentile]; m->mothurOut("Optimizing maxhomop to " + toString(maxHomoP) + "."); m->mothurOutEndLine(); }
			else if (optimize[i] == "minlength") { int mincriteriaPercentile = int(seqLength.size() * ((100 - criteria) / (float) 100)); minLength = seqLength[mincriteriaPercentile]; m->mothurOut("Optimizing minlength to " + toString(minLength) + "."); m->mothurOutEndLine(); }
			else if (optimize[i] == "maxlength") { maxLength = seqLength[criteriaPercentile]; m->mothurOut("Optimizing maxlength to " + toString(maxLength) + "."); m->mothurOutEndLine(); }
		}
        
#ifdef USE_MPI
    }
    
    MPI_Status status; 
    MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
    MPI_Comm_size(MPI_COMM_WORLD, &processors); 
    
    if (pid == 0) { 
        //send file positions to all processes
        for(int i = 1; i < processors; i++) { 
            MPI_Send(&startPos, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
            MPI_Send(&endPos, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
            MPI_Send(&maxAmbig, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
            MPI_Send(&maxHomoP, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
            MPI_Send(&minLength, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
            MPI_Send(&maxLength, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
        }
    }else {
        MPI_Recv(&startPos, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
        MPI_Recv(&endPos, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
        MPI_Recv(&maxAmbig, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
        MPI_Recv(&maxHomoP, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
        MPI_Recv(&minLength, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
        MPI_Recv(&maxLength, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
    }
    MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
#endif
        return 0;
        
    }
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "getSummaryReport");
		exit(1);
	}
}
//***************************************************************************************************************
int ScreenSeqsCommand::optimizeContigs(){
	try {
		vector<int> olengths;
		vector<int> oStarts;
		vector<int> oEnds;
		vector<int> numMismatches;
        vector<int> numNs;
		
        vector<unsigned long long> positions;
        vector<linePair> contigsLines;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		positions = m->divideFilePerLine(contigsreport, processors);
		for (int i = 0; i < (positions.size()-1); i++) { contigsLines.push_back(linePair(positions[i], positions[(i+1)])); }	
#else
		if(processors == 1){ contigsLines.push_back(linePair(0, 1000));  }
        else {
            int numContigsSeqs = 0;
            positions = m->setFilePosEachLine(contigsreport, numContigsSeqs); 
            if (positions.size() < processors) { processors = positions.size(); }
            
            //figure out how many sequences you have to process
            int numSeqsPerProcessor = numContigsSeqs / processors;
            for (int i = 0; i < processors; i++) {
                int startIndex =  i * numSeqsPerProcessor;
                if(i == (processors - 1)){	numSeqsPerProcessor = numContigsSeqs - i * numSeqsPerProcessor; 	}
                contigsLines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
            }
        }
#endif
		
#ifdef USE_MPI
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
		
		if (pid == 0) { 
			driverContigsSummary(olengths, oStarts, oEnds, numMismatches, numNs, contigsLines[0]);
#else
            createProcessesContigsSummary(olengths, oStarts, oEnds, numMismatches, numNs, contigsLines); 
            
			if (m->control_pressed) {  return 0; }
#endif
            sort(olengths.begin(), olengths.end());
            sort(oStarts.begin(), oStarts.end());
            sort(oEnds.begin(), oEnds.end());
            sort(numMismatches.begin(), numMismatches.end());
            sort(numNs.begin(), numNs.end());
            
            //numSeqs is the number of unique seqs, startPosition.size() is the total number of seqs, we want to optimize using all seqs
            int criteriaPercentile	= int(oStarts.size() * (criteria / (float) 100));
            
            for (int i = 0; i < optimize.size(); i++) {
                if (optimize[i] == "ostart") { oStart = oStarts[criteriaPercentile]; m->mothurOut("Optimizing ostart to " + toString(oStart) + "."); m->mothurOutEndLine(); }
                else if (optimize[i] == "oend") { int endcriteriaPercentile = int(oEnds.size() * ((100 - criteria) / (float) 100));  oEnd = oEnds[endcriteriaPercentile]; m->mothurOut("Optimizing oend to " + toString(oEnd) + "."); m->mothurOutEndLine();}
                else if (optimize[i] == "mismatches") { mismatches = numMismatches[criteriaPercentile]; m->mothurOut("Optimizing mismatches to " + toString(mismatches) + "."); m->mothurOutEndLine(); }
                else if (optimize[i] == "maxn") { maxN = numNs[criteriaPercentile]; m->mothurOut("Optimizing maxn to " + toString(maxN) + "."); m->mothurOutEndLine(); }
                else if (optimize[i] == "minoverlap") { int mincriteriaPercentile = int(olengths.size() * ((100 - criteria) / (float) 100)); minOverlap = olengths[mincriteriaPercentile]; m->mothurOut("Optimizing minoverlap to " + toString(minOverlap) + "."); m->mothurOutEndLine(); }

            }
            
#ifdef USE_MPI
		}
		
		MPI_Status status; 
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
		MPI_Comm_size(MPI_COMM_WORLD, &processors); 
        
		if (pid == 0) { 
			//send file positions to all processes
			for(int i = 1; i < processors; i++) { 
                MPI_Send(&minOverlap, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
				MPI_Send(&oStart, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
				MPI_Send(&oEnd, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
				MPI_Send(&mismatches, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
                MPI_Send(&maxN, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
			}
		}else {
            MPI_Recv(&minOverlap, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
			MPI_Recv(&oStart, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
			MPI_Recv(&oEnd, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
			MPI_Recv(&mismatches, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
            MPI_Recv(&maxN, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
		}
		MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
#endif
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "optimizeContigs");
		exit(1);
	}
}
/**************************************************************************************/
int ScreenSeqsCommand::driverContigsSummary(vector<int>& oLength, vector<int>& ostartPosition, vector<int>& oendPosition, vector<int>& omismatches, vector<int>& numNs, linePair filePos) {	
	try {
		
        string name;
        //Name	Length	Overlap_Length	Overlap_Start	Overlap_End	MisMatches	Num_Ns
        int length, OLength, thisOStart, thisOEnd, numMisMatches, numns;
        
  		ifstream in;
		m->openInputFile(contigsreport, in);
        
		in.seekg(filePos.start);
        if (filePos.start == 0) { //read headers
            m->zapGremlins(in); m->gobble(in); m->getline(in); m->gobble(in);
        }
        
		bool done = false;
		int count = 0;
        
		while (!done) {
            
			if (m->control_pressed) { in.close(); return 1; }
            
            //seqname	start	end	nbases	ambigs	polymer	numSeqs
            in >> name >> length >> OLength >> thisOStart >> thisOEnd >> numMisMatches >> numns; m->gobble(in);
            
            int num = 1;
            if ((namefile != "") || (countfile !="")){
                //make sure this sequence is in the namefile, else error 
                map<string, int>::iterator it = nameMap.find(name);
                
                if (it == nameMap.end()) { m->mothurOut("[ERROR]: " + name + " is not in your namefile, please correct."); m->mothurOutEndLine(); m->control_pressed = true; }
                else { num = it->second; }
            }
            
            //for each sequence this sequence represents
            for (int i = 0; i < num; i++) {
                ostartPosition.push_back(thisOStart);
                oendPosition.push_back(thisOEnd);
                oLength.push_back(OLength);
                omismatches.push_back(numMisMatches);
                numNs.push_back(numns);
            }
            
            count++;
			
			//if((count) % 100 == 0){	m->mothurOut("Optimizing sequence: " + toString(count)); m->mothurOutEndLine();		}
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
            unsigned long long pos = in.tellg();
            if ((pos == -1) || (pos >= filePos.end)) { break; }
#else
            if (in.eof()) { break; }
#endif
		}
		
		in.close();
		
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "driverContigsSummary");
		exit(1);
	}
}

/**************************************************************************************************/
int ScreenSeqsCommand::createProcessesContigsSummary(vector<int>& oLength, vector<int>& ostartPosition, vector<int>& oendPosition, vector<int>& omismatches, vector<int>& numNs, vector<linePair> contigsLines) {
	try {
        
        int process = 1;
		int num = 0;
		vector<int> processIDS;
        bool recalc = false;
        
        
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driverContigsSummary(oLength, ostartPosition, oendPosition, omismatches, numNs, contigsLines[process]);
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = contigsreport + m->mothurGetpid(process) + ".num.temp";
				m->openOutputFile(tempFile, out);
				
				out << num << endl;
				out << ostartPosition.size() << endl;
				for (int k = 0; k < ostartPosition.size(); k++)		{		out << ostartPosition[k] << '\t';   }  out << endl;
				for (int k = 0; k < oendPosition.size(); k++)		{		out << oendPosition[k] << '\t';     }  out << endl;
				for (int k = 0; k < oLength.size(); k++)			{		out << oLength[k] << '\t';          }  out << endl;
				for (int k = 0; k < omismatches.size(); k++)        {		out << omismatches[k] << '\t';      }  out << endl;
                for (int k = 0; k < numNs.size(); k++)              {		out << numNs[k] << '\t';            }  out << endl;
				
				out.close();
				
				exit(0);
			}else { 
                m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(process) + "\n"); processors = process;
                for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                //wait to die
                for (int i=0;i<processIDS.size();i++) {
                    int temp = processIDS[i];
                    wait(&temp);
                }
                m->control_pressed = false;
                for (int i=0;i<processIDS.size();i++) {
                    m->mothurRemove(contigsreport + (toString(processIDS[i]) + ".num.temp"));
                }
                recalc = true;
                break;
			}
		}
        
        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->control_pressed = false;  for (int i=0;i<processIDS.size();i++) {m->mothurRemove(contigsreport + (toString(processIDS[i]) + ".num.temp"));}processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
            
            vector<unsigned long long> positions;
            contigsLines.clear();
            positions = m->divideFilePerLine(contigsreport, processors);
            for (int i = 0; i < (positions.size()-1); i++) { contigsLines.push_back(linePair(positions[i], positions[(i+1)])); }
            
            //redo file divide
            lines.clear(); positions.clear();
            positions = m->divideFile(fastafile, processors);
            for (int i = 0; i < (positions.size()-1); i++) {  lines.push_back(linePair(positions[i], positions[(i+1)]));  }
            
            ostartPosition.clear();
            oendPosition.clear();
            oLength.clear();
            omismatches.clear();
            numNs.clear();
            
            num = 0;
            processIDS.resize(0);
            process = 1;
            
            //loop through and create all the processes you want
            while (process != processors) {
                pid_t pid = fork();
                
                if (pid > 0) {
                    processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                    process++;
                }else if (pid == 0){
                    num = driverContigsSummary(oLength, ostartPosition, oendPosition, omismatches, numNs, contigsLines[process]);
                    
                    //pass numSeqs to parent
                    ofstream out;
                    string tempFile = contigsreport + m->mothurGetpid(process) + ".num.temp";
                    m->openOutputFile(tempFile, out);
                    
                    out << num << endl;
                    out << ostartPosition.size() << endl;
                    for (int k = 0; k < ostartPosition.size(); k++)		{		out << ostartPosition[k] << '\t';   }  out << endl;
                    for (int k = 0; k < oendPosition.size(); k++)		{		out << oendPosition[k] << '\t';     }  out << endl;
                    for (int k = 0; k < oLength.size(); k++)			{		out << oLength[k] << '\t';          }  out << endl;
                    for (int k = 0; k < omismatches.size(); k++)        {		out << omismatches[k] << '\t';      }  out << endl;
                    for (int k = 0; k < numNs.size(); k++)              {		out << numNs[k] << '\t';            }  out << endl;
                    
                    out.close();
                    
                    exit(0);
                }else { 
                    m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
                    for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                    exit(0);
                }
            }
        }
		
		num = driverContigsSummary(oLength, ostartPosition, oendPosition, omismatches, numNs, contigsLines[0]);
        
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		//parent reads in and combine Filter info
		for (int i = 0; i < processIDS.size(); i++) {
			string tempFilename = contigsreport + toString(processIDS[i]) + ".num.temp";
			ifstream in;
			m->openInputFile(tempFilename, in);
			
			int temp, tempNum;
			in >> tempNum; m->gobble(in); num += tempNum;
			in >> tempNum; m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{		in >> temp; ostartPosition.push_back(temp);		}		m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{		in >> temp; oendPosition.push_back(temp);		}		m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{		in >> temp; oLength.push_back(temp);			}		m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{		in >> temp; omismatches.push_back(temp);        }		m->gobble(in);
            for (int k = 0; k < tempNum; k++)			{		in >> temp; numNs.push_back(temp);              }		m->gobble(in);
            
			in.close();
			m->mothurRemove(tempFilename);
		}
		
		
#else 
        //////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the seqSumData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//Taking advantage of shared memory to allow both threads to add info to vectors.
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		/*
		vector<contigsSumData*> pDataArray; 
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1]; 
		
		//Create processor worker threads.
		for( int i=0; i<processors-1; i++ ){
            
			// Allocate memory for thread data.
			contigsSumData* tempSum = new contigsSumData(contigsreport, m, contigsLines[i].start, contigsLines[i].end, namefile, countfile, nameMap);
			pDataArray.push_back(tempSum);
			
			//MySeqSumThreadFunction is in header. It must be global or static to work with the threads.
			//default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
			hThreadArray[i] = CreateThread(NULL, 0, MyContigsSumThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);   
		}
		*/
        contigsLines[processors-1].start = 0;
        //do your part
		num = driverContigsSummary(oLength, ostartPosition, oendPosition, omismatches, numNs, contigsLines[processors-1]);
        /*
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
			num += pDataArray[i]->count;
            if (pDataArray[i]->count != pDataArray[i]->end) {
                m->mothurOut("[ERROR]: process " + toString(i) + " only processed " + toString(pDataArray[i]->count) + " of " + toString(pDataArray[i]->end) + " sequences assigned to it, quitting. \n"); m->control_pressed = true; 
            }
            for (int k = 0; k < pDataArray[i]->ostartPosition.size(); k++)  {	ostartPosition.push_back(pDataArray[i]->ostartPosition[k]);     }
			for (int k = 0; k < pDataArray[i]->oendPosition.size(); k++)    {	oendPosition.push_back(pDataArray[i]->oendPosition[k]);         }
            for (int k = 0; k < pDataArray[i]->oLength.size(); k++)         {	oLength.push_back(pDataArray[i]->oLength[k]);                   }
            for (int k = 0; k < pDataArray[i]->omismatches.size(); k++)     {	omismatches.push_back(pDataArray[i]->omismatches[k]);           }
            for (int k = 0; k < pDataArray[i]->numNs.size(); k++)           {	numNs.push_back(pDataArray[i]->numNs[k]);                       }
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
        */
#endif		
        return num;
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "createProcessesContigsSummary");
		exit(1);
	}
}
//***************************************************************************************************************
int ScreenSeqsCommand::optimizeAlign(){
	try {
        
		vector<float> sims;
		vector<float> scores;
		vector<int> inserts;
		
        vector<unsigned long long> positions;
        vector<linePair> alignLines;
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		positions = m->divideFilePerLine(alignreport, processors);
		for (int i = 0; i < (positions.size()-1); i++) { alignLines.push_back(linePair(positions[i], positions[(i+1)])); }	
#else
		if(processors == 1){ alignLines.push_back(linePair(0, 1000));  }
        else {
            int numAlignSeqs = 0;
            positions = m->setFilePosEachLine(alignreport, numAlignSeqs); 
            if (positions.size() < processors) { processors = positions.size(); }
            
            //figure out how many sequences you have to process
            int numSeqsPerProcessor = numAlignSeqs / processors;
            for (int i = 0; i < processors; i++) {
                int startIndex =  i * numSeqsPerProcessor;
                if(i == (processors - 1)){	numSeqsPerProcessor = numAlignSeqs - i * numSeqsPerProcessor; 	}
                alignLines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
            }
        }
#endif
		
#ifdef USE_MPI
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
		
		if (pid == 0) { 
			driverAlignSummary(sims, scores, inserts, alignLines[0]);
#else
            createProcessesAlignSummary(sims, scores, inserts, alignLines); 
            
			if (m->control_pressed) {  return 0; }
#endif
            sort(sims.begin(), sims.end());
            sort(scores.begin(), scores.end());
            sort(inserts.begin(), inserts.end());
            
            //numSeqs is the number of unique seqs, startPosition.size() is the total number of seqs, we want to optimize using all seqs
            int criteriaPercentile	= int(sims.size() * (criteria / (float) 100));
            
            for (int i = 0; i < optimize.size(); i++) {
                if (optimize[i] == "minsim") { int mincriteriaPercentile = int(sims.size() * ((100 - criteria) / (float) 100)); minSim = sims[mincriteriaPercentile];  m->mothurOut("Optimizing minsim to " + toString(minSim) + "."); m->mothurOutEndLine();}
                else if (optimize[i] == "minscore") { int mincriteriaPercentile = int(scores.size() * ((100 - criteria) / (float) 100)); minScore = scores[mincriteriaPercentile];  m->mothurOut("Optimizing minscore to " + toString(minScore) + "."); m->mothurOutEndLine(); }
                else if (optimize[i] == "maxinsert") { maxInsert = inserts[criteriaPercentile]; m->mothurOut("Optimizing maxinsert to " + toString(maxInsert) + "."); m->mothurOutEndLine(); }
            }
            
#ifdef USE_MPI
		}
		
		MPI_Status status; 
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
		MPI_Comm_size(MPI_COMM_WORLD, &processors); 
        
		if (pid == 0) { 
			//send file positions to all processes
			for(int i = 1; i < processors; i++) { 
                MPI_Send(&minSim, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
				MPI_Send(&minScore, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
				MPI_Send(&maxInsert, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
			}
		}else {
            MPI_Recv(&minSim, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
			MPI_Recv(&minScore, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
			MPI_Recv(&maxInsert, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
		}
		MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
#endif
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "optimizeContigs");
		exit(1);
	}
}
/**************************************************************************************/
int ScreenSeqsCommand::driverAlignSummary(vector<float>& sims, vector<float>& scores, vector<int>& inserts, linePair filePos) {	
	try {
		
        string name, TemplateName, SearchMethod, AlignmentMethod;
        //QueryName	QueryLength	TemplateName	TemplateLength	SearchMethod	SearchScore	AlignmentMethod	QueryStart	QueryEnd	TemplateStart	TemplateEnd	PairwiseAlignmentLength	GapsInQuery	GapsInTemplate	LongestInsert	SimBtwnQuery&Template
        //checking for minScore, maxInsert, minSim
        int length, TemplateLength,	 QueryStart,	QueryEnd,	TemplateStart,	TemplateEnd,	PairwiseAlignmentLength,	GapsInQuery,	GapsInTemplate,	LongestInsert;
        float SearchScore, SimBtwnQueryTemplate;
         
  		ifstream in;
		m->openInputFile(alignreport, in);
        
		in.seekg(filePos.start);
        if (filePos.start == 0) { //read headers
            m->zapGremlins(in); m->gobble(in); m->getline(in); m->gobble(in);
        }
        
		bool done = false;
		int count = 0;
        
		while (!done) {
            
			if (m->control_pressed) { in.close(); return 1; }
            
            in >> name >> length >> TemplateName >> TemplateLength >> SearchMethod >> SearchScore >> AlignmentMethod >> QueryStart >> QueryEnd >> TemplateStart >> TemplateEnd >> PairwiseAlignmentLength >> GapsInQuery >> GapsInTemplate >> LongestInsert >> SimBtwnQueryTemplate; m->gobble(in);
            
            int num = 1;
            if ((namefile != "") || (countfile !="")){
                //make sure this sequence is in the namefile, else error 
                map<string, int>::iterator it = nameMap.find(name);
                
                if (it == nameMap.end()) { m->mothurOut("[ERROR]: " + name + " is not in your namefile, please correct."); m->mothurOutEndLine(); m->control_pressed = true; }
                else { num = it->second; }
            }
            
            //for each sequence this sequence represents
            for (int i = 0; i < num; i++) {
                sims.push_back(SimBtwnQueryTemplate);
                scores.push_back(SearchScore);
                inserts.push_back(LongestInsert);
            }
            
            count++;
			
			//if((count) % 100 == 0){	m->mothurOut("Optimizing sequence: " + toString(count)); m->mothurOutEndLine();		}
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
            unsigned long long pos = in.tellg();
            if ((pos == -1) || (pos >= filePos.end)) { break; }
#else
            if (in.eof()) { break; }
#endif
		}
		
		in.close();
		
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "driverAlignSummary");
		exit(1);
	}
}

/**************************************************************************************************/
int ScreenSeqsCommand::createProcessesAlignSummary(vector<float>& sims, vector<float>& scores, vector<int>& inserts, vector<linePair> alignLines) {
	try {
        
        int process = 1;
		int num = 0;
		vector<int> processIDS;
        bool recalc = false;
        
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
        
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driverAlignSummary(sims, scores, inserts, alignLines[process]);
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = alignreport + m->mothurGetpid(process) + ".num.temp";
				m->openOutputFile(tempFile, out);
				
				out << num << endl;
				out << sims.size() << endl;
				for (int k = 0; k < sims.size(); k++)		{		out << sims[k] << '\t';         }  out << endl;
				for (int k = 0; k < scores.size(); k++)		{		out << scores[k] << '\t';       }  out << endl;
				for (int k = 0; k < inserts.size(); k++)	{		out << inserts[k] << '\t';      }  out << endl;
				
				out.close();
				
				exit(0);
			}else { 
                m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(process) + "\n"); processors = process;
                for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                //wait to die
                for (int i=0;i<processIDS.size();i++) {
                    int temp = processIDS[i];
                    wait(&temp);
                }
                m->control_pressed = false;
                for (int i=0;i<processIDS.size();i++) {
                    m->mothurRemove(contigsreport + (toString(processIDS[i]) + ".num.temp"));
                }
                recalc = true;
                break;
			}
		}
        
        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->control_pressed = false;  for (int i=0;i<processIDS.size();i++) {m->mothurRemove(contigsreport + (toString(processIDS[i]) + ".num.temp"));}processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
            
            vector<unsigned long long> positions;
            alignLines.clear();
            positions = m->divideFilePerLine(contigsreport, processors);
            for (int i = 0; i < (positions.size()-1); i++) { alignLines.push_back(linePair(positions[i], positions[(i+1)])); }
            
            //redo file divide
            lines.clear(); positions.clear();
            positions = m->divideFile(fastafile, processors);
            for (int i = 0; i < (positions.size()-1); i++) {  lines.push_back(linePair(positions[i], positions[(i+1)]));  }
            
            sims.clear();
            scores.clear();
            inserts.clear();
            
            num = 0;
            processIDS.resize(0);
            process = 1;
            
            //loop through and create all the processes you want
            while (process != processors) {
                pid_t pid = fork();
                
                if (pid > 0) {
                    processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                    process++;
                }else if (pid == 0){
                    num = driverAlignSummary(sims, scores, inserts, alignLines[process]);
                    
                    //pass numSeqs to parent
                    ofstream out;
                    string tempFile = alignreport + m->mothurGetpid(process) + ".num.temp";
                    m->openOutputFile(tempFile, out);
                    
                    out << num << endl;
                    out << sims.size() << endl;
                    for (int k = 0; k < sims.size(); k++)		{		out << sims[k] << '\t';         }  out << endl;
                    for (int k = 0; k < scores.size(); k++)		{		out << scores[k] << '\t';       }  out << endl;
                    for (int k = 0; k < inserts.size(); k++)	{		out << inserts[k] << '\t';      }  out << endl;
                    
                    out.close();
                    
                    exit(0);
                }else { 
                    m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
                    for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                    exit(0);
                }
            }
        }
		
		num = driverAlignSummary(sims, scores, inserts, alignLines[0]);
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		//parent reads in and combine Filter info
		for (int i = 0; i < processIDS.size(); i++) {
			string tempFilename = alignreport + toString(processIDS[i]) + ".num.temp";
			ifstream in;
			m->openInputFile(tempFilename, in);
			
			int temp, tempNum;
            float temp2;
			in >> tempNum; m->gobble(in); num += tempNum;
			in >> tempNum; m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{		in >> temp2; sims.push_back(temp2);		}		m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{		in >> temp2; scores.push_back(temp2);		}		m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{		in >> temp; inserts.push_back(temp);	}		m->gobble(in);
			  
			in.close();
			m->mothurRemove(tempFilename);
		}
		
		
#else 
        //////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the seqSumData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//Taking advantage of shared memory to allow both threads to add info to vectors.
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		/*
		vector<alignsData*> pDataArray; 
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1]; 
		
		//Create processor worker threads.
		for( int i=0; i<processors-1; i++ ){
            
			// Allocate memory for thread data.
			alignsData* tempSum = new alignsData(alignreport, m, alignLines[i].start, alignLines[i].end, namefile, countfile, nameMap);
			pDataArray.push_back(tempSum);
			
			//MySeqSumThreadFunction is in header. It must be global or static to work with the threads.
			//default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
			hThreadArray[i] = CreateThread(NULL, 0, MyAlignsThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);   
		}*/
		alignLines[processors-1].start = 0;
        //do your part
		num = driverAlignSummary(sims, scores, inserts, alignLines[processors-1]);
       /*
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
			num += pDataArray[i]->count;
            if (pDataArray[i]->count != pDataArray[i]->end) {
                m->mothurOut("[ERROR]: process " + toString(i) + " only processed " + toString(pDataArray[i]->count) + " of " + toString(pDataArray[i]->end) + " sequences assigned to it, quitting. \n"); m->control_pressed = true; 
            }
            for (int k = 0; k < pDataArray[i]->sims.size(); k++)        {	sims.push_back(pDataArray[i]->sims[k]);         }
			for (int k = 0; k < pDataArray[i]->scores.size(); k++)      {	scores.push_back(pDataArray[i]->scores[k]);     }
            for (int k = 0; k < pDataArray[i]->inserts.size(); k++)     {	inserts.push_back(pDataArray[i]->inserts[k]);   }
           	CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
        */
#endif		
        return num;
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "createProcessesAlignSummary");
		exit(1);
	}
}
//***************************************************************************************************************
int ScreenSeqsCommand::getSummary(vector<unsigned long long>& positions){
	try {
		
		vector<int> startPosition;
		vector<int> endPosition;
		vector<int> seqLength;
		vector<int> ambigBases;
		vector<int> longHomoPolymer;
        vector<int> numNs;
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
		positions = m->divideFile(fastafile, processors);
		for (int i = 0; i < (positions.size()-1); i++) { lines.push_back(linePair(positions[i], positions[(i+1)])); }	
#else
		if(processors == 1){ lines.push_back(linePair(0, 1000));  }
        else {
            int numFastaSeqs = 0;
            positions = m->setFilePosFasta(fastafile, numFastaSeqs); 
            if (positions.size() < processors) { processors = positions.size(); }
            
            //figure out how many sequences you have to process
            int numSeqsPerProcessor = numFastaSeqs / processors;
            for (int i = 0; i < processors; i++) {
                int startIndex =  i * numSeqsPerProcessor;
                if(i == (processors - 1)){	numSeqsPerProcessor = numFastaSeqs - i * numSeqsPerProcessor; 	}
                lines.push_back(linePair(positions[startIndex], numSeqsPerProcessor));
            }
        }
#endif
		
#ifdef USE_MPI
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
		
		if (pid == 0) {
            linePair tempLine(0, positions[positions.size()-1]);
			driverCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, numNs, fastafile, tempLine);
#else
		int numSeqs = 0;
		//#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
			if(processors == 1){
				numSeqs = driverCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, numNs, fastafile, lines[0]);
			}else{
				numSeqs = createProcessesCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, numNs, fastafile); 
			}
				
			if (m->control_pressed) {  return 0; }
#endif
		sort(startPosition.begin(), startPosition.end());
		sort(endPosition.begin(), endPosition.end());
		sort(seqLength.begin(), seqLength.end());
		sort(ambigBases.begin(), ambigBases.end());
		sort(longHomoPolymer.begin(), longHomoPolymer.end());
        sort(numNs.begin(), numNs.end());
		
		//numSeqs is the number of unique seqs, startPosition.size() is the total number of seqs, we want to optimize using all seqs
		int criteriaPercentile	= int(startPosition.size() * (criteria / (float) 100));
		
		for (int i = 0; i < optimize.size(); i++) {
			if (optimize[i] == "start") { startPos = startPosition[criteriaPercentile]; m->mothurOut("Optimizing start to " + toString(startPos) + "."); m->mothurOutEndLine(); }
			else if (optimize[i] == "end") { int endcriteriaPercentile = int(endPosition.size() * ((100 - criteria) / (float) 100));  endPos = endPosition[endcriteriaPercentile]; m->mothurOut("Optimizing end to " + toString(endPos) + "."); m->mothurOutEndLine();}
			else if (optimize[i] == "maxambig") { maxAmbig = ambigBases[criteriaPercentile]; m->mothurOut("Optimizing maxambig to " + toString(maxAmbig) + "."); m->mothurOutEndLine(); }
			else if (optimize[i] == "maxhomop") { maxHomoP = longHomoPolymer[criteriaPercentile]; m->mothurOut("Optimizing maxhomop to " + toString(maxHomoP) + "."); m->mothurOutEndLine(); }
			else if (optimize[i] == "minlength") { int mincriteriaPercentile = int(seqLength.size() * ((100 - criteria) / (float) 100)); minLength = seqLength[mincriteriaPercentile]; m->mothurOut("Optimizing minlength to " + toString(minLength) + "."); m->mothurOutEndLine(); }
			else if (optimize[i] == "maxlength") { maxLength = seqLength[criteriaPercentile]; m->mothurOut("Optimizing maxlength to " + toString(maxLength) + "."); m->mothurOutEndLine(); }
            else if (optimize[i] == "maxn") { maxN = numNs[criteriaPercentile]; m->mothurOut("Optimizing maxn to " + toString(maxN) + "."); m->mothurOutEndLine(); }
		}

#ifdef USE_MPI
		}
		
		MPI_Status status; 
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); 
		MPI_Comm_size(MPI_COMM_WORLD, &processors); 
			
		if (pid == 0) { 
			//send file positions to all processes
			for(int i = 1; i < processors; i++) { 
				MPI_Send(&startPos, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
				MPI_Send(&endPos, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
				MPI_Send(&maxAmbig, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
				MPI_Send(&maxHomoP, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
				MPI_Send(&minLength, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
				MPI_Send(&maxLength, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
                MPI_Send(&maxN, 1, MPI_INT, i, 2001, MPI_COMM_WORLD);
			}
		}else {
			MPI_Recv(&startPos, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
			MPI_Recv(&endPos, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
			MPI_Recv(&maxAmbig, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
			MPI_Recv(&maxHomoP, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
			MPI_Recv(&minLength, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
			MPI_Recv(&maxLength, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
            MPI_Recv(&maxN, 1, MPI_INT, 0, 2001, MPI_COMM_WORLD, &status);
		}
		MPI_Barrier(MPI_COMM_WORLD); //make everyone wait - just in case
#endif
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "getSummary");
		exit(1);
	}
}
/**************************************************************************************/
int ScreenSeqsCommand::driverCreateSummary(vector<int>& startPosition, vector<int>& endPosition, vector<int>& seqLength, vector<int>& ambigBases, vector<int>& longHomoPolymer, vector<int>& numNs, string filename, linePair filePos) {	
	try {
		
		ifstream in;
		m->openInputFile(filename, in);
				
		in.seekg(filePos.start);
        
        //adjust start if null strings
        if (filePos.start == 0) {  m->zapGremlins(in); m->gobble(in);  }
        
        m->gobble(in);

		bool done = false;
		int count = 0;
	
		while (!done) {
				
			if (m->control_pressed) { in.close(); return 1; }
					
			Sequence current(in); m->gobble(in);
	
			if (current.getName() != "") {
				int num = 1;
				if ((namefile != "") || (countfile !="")){
					//make sure this sequence is in the namefile, else error 
					map<string, int>::iterator it = nameMap.find(current.getName());
					
					if (it == nameMap.end()) { m->mothurOut("[ERROR]: " + current.getName() + " is not in your namefile, please correct."); m->mothurOutEndLine(); m->control_pressed = true; }
					else { num = it->second; }
				}
				
				//for each sequence this sequence represents
                int numns = current.getNumNs();
				for (int i = 0; i < num; i++) {
					startPosition.push_back(current.getStartPos());
					endPosition.push_back(current.getEndPos());
					seqLength.push_back(current.getNumBases());
					ambigBases.push_back(current.getAmbigBases());
					longHomoPolymer.push_back(current.getLongHomoPolymer());
                    numNs.push_back(numns);
				}
				
				count++;
			}
			//if((count) % 100 == 0){	m->mothurOut("Optimizing sequence: " + toString(count)); m->mothurOutEndLine();		}
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				unsigned long long pos = in.tellg();
				if ((pos == -1) || (pos >= filePos.end)) { break; }
			#else
				if (in.eof()) { break; }
			#endif
			
		}
		
		in.close();
		
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "driverCreateSummary");
		exit(1);
	}
}
/**************************************************************************************************/
int ScreenSeqsCommand::createProcessesCreateSummary(vector<int>& startPosition, vector<int>& endPosition, vector<int>& seqLength, vector<int>& ambigBases, vector<int>& longHomoPolymer, vector<int>& numNs, string filename) {
	try {
        
        int process = 1;
		int num = 0;
		vector<int> processIDS;
        bool recalc = false;

#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driverCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, numNs, fastafile, lines[process]);
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = fastafile + m->mothurGetpid(process) + ".num.temp";
				m->openOutputFile(tempFile, out);
				
				out << num << endl;
				out << startPosition.size() << endl;
				for (int k = 0; k < startPosition.size(); k++)		{		out << startPosition[k] << '\t'; }  out << endl;
				for (int k = 0; k < endPosition.size(); k++)		{		out << endPosition[k] << '\t'; }  out << endl;
				for (int k = 0; k < seqLength.size(); k++)			{		out << seqLength[k] << '\t'; }  out << endl;
				for (int k = 0; k < ambigBases.size(); k++)			{		out << ambigBases[k] << '\t'; }  out << endl;
				for (int k = 0; k < longHomoPolymer.size(); k++)	{		out << longHomoPolymer[k] << '\t'; }  out << endl;
                for (int k = 0; k < numNs.size(); k++)	{		out << numNs[k] << '\t'; }  out << endl;
				
				out.close();
				
				exit(0);
			}else { 
                m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(process) + "\n"); processors = process;
                for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                //wait to die
                for (int i=0;i<processIDS.size();i++) {
                    int temp = processIDS[i];
                    wait(&temp);
                }
                m->control_pressed = false;
                for (int i=0;i<processIDS.size();i++) {
                    m->mothurRemove(contigsreport + (toString(processIDS[i]) + ".num.temp"));
                }
                recalc = true;
                break;
			}
		}
        
        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->control_pressed = false;  for (int i=0;i<processIDS.size();i++) {m->mothurRemove(contigsreport + (toString(processIDS[i]) + ".num.temp"));}processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
            
            vector<unsigned long long> positions;
            lines.clear();
            positions = m->divideFile(fastafile, processors);
            for (int i = 0; i < (positions.size()-1); i++) {  lines.push_back(linePair(positions[i], positions[(i+1)]));  }
            
            startPosition.clear();
            endPosition.clear();
            seqLength.clear();
            ambigBases.clear();
            longHomoPolymer.clear();
            numNs.clear();
            
            num = 0;
            processIDS.resize(0);
            process = 1;
            
            //loop through and create all the processes you want
            while (process != processors) {
                pid_t pid = fork();
                
                if (pid > 0) {
                    processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                    process++;
                }else if (pid == 0){
                    num = driverCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, numNs, fastafile, lines[process]);
                    
                    //pass numSeqs to parent
                    ofstream out;
                    string tempFile = fastafile + m->mothurGetpid(process) + ".num.temp";
                    m->openOutputFile(tempFile, out);
                    
                    out << num << endl;
                    out << startPosition.size() << endl;
                    for (int k = 0; k < startPosition.size(); k++)		{		out << startPosition[k] << '\t'; }  out << endl;
                    for (int k = 0; k < endPosition.size(); k++)		{		out << endPosition[k] << '\t'; }  out << endl;
                    for (int k = 0; k < seqLength.size(); k++)			{		out << seqLength[k] << '\t'; }  out << endl;
                    for (int k = 0; k < ambigBases.size(); k++)			{		out << ambigBases[k] << '\t'; }  out << endl;
                    for (int k = 0; k < longHomoPolymer.size(); k++)	{		out << longHomoPolymer[k] << '\t'; }  out << endl;
                    for (int k = 0; k < numNs.size(); k++)	{		out << numNs[k] << '\t'; }  out << endl;
                    
                    out.close();
                    
                    exit(0);
                }else { 
                    m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
                    for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                    exit(0);
                }
            }
        }

		
		num = driverCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, numNs, fastafile, lines[0]);
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		//parent reads in and combine Filter info
		for (int i = 0; i < processIDS.size(); i++) {
			string tempFilename = fastafile + toString(processIDS[i]) + ".num.temp";
			ifstream in;
			m->openInputFile(tempFilename, in);
			
			int temp, tempNum;
			in >> tempNum; m->gobble(in); num += tempNum;
			in >> tempNum; m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{		in >> temp; startPosition.push_back(temp);		}		m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{		in >> temp; endPosition.push_back(temp);		}		m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{		in >> temp; seqLength.push_back(temp);			}		m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{		in >> temp; ambigBases.push_back(temp);			}		m->gobble(in);
			for (int k = 0; k < tempNum; k++)			{		in >> temp; longHomoPolymer.push_back(temp);	}		m->gobble(in);
            for (int k = 0; k < tempNum; k++)			{		in >> temp; numNs.push_back(temp);	}		m->gobble(in);
				
			in.close();
			m->mothurRemove(tempFilename);
		}
		
		
#else 
        //////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the seqSumData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//Taking advantage of shared memory to allow both threads to add info to vectors.
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		vector<sumData*> pDataArray; 
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1]; 
		
		//Create processor worker threads.
		for( int i=0; i<processors-1; i++ ){
            
			// Allocate memory for thread data.
			sumData* tempSum = new sumData(filename, m, lines[i].start, lines[i].end, namefile, countfile, nameMap);
			pDataArray.push_back(tempSum);
			
			//MySeqSumThreadFunction is in header. It must be global or static to work with the threads.
			//default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
			hThreadArray[i] = CreateThread(NULL, 0, MySumThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);   
		}
		
        //do your part
		num = driverCreateSummary(startPosition, endPosition, seqLength, ambigBases, longHomoPolymer, numNs, fastafile, lines[processors-1]);
         
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
			num += pDataArray[i]->count;
            if (pDataArray[i]->count != pDataArray[i]->end) {
                m->mothurOut("[ERROR]: process " + toString(i) + " only processed " + toString(pDataArray[i]->count) + " of " + toString(pDataArray[i]->end) + " sequences assigned to it, quitting. \n"); m->control_pressed = true; 
            }
            for (int k = 0; k < pDataArray[i]->startPosition.size(); k++) {	startPosition.push_back(pDataArray[i]->startPosition[k]);       }
			for (int k = 0; k < pDataArray[i]->endPosition.size(); k++) {	endPosition.push_back(pDataArray[i]->endPosition[k]);       }
            for (int k = 0; k < pDataArray[i]->seqLength.size(); k++) {	seqLength.push_back(pDataArray[i]->seqLength[k]);       }
            for (int k = 0; k < pDataArray[i]->ambigBases.size(); k++) {	ambigBases.push_back(pDataArray[i]->ambigBases[k]);       }
            for (int k = 0; k < pDataArray[i]->longHomoPolymer.size(); k++) {	longHomoPolymer.push_back(pDataArray[i]->longHomoPolymer[k]);       }
            for (int k = 0; k < pDataArray[i]->numNs.size(); k++) {	numNs.push_back(pDataArray[i]->numNs[k]);       }
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}

#endif		
        return num;
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "createProcessesCreateSummary");
		exit(1);
	}
}

//***************************************************************************************************************

int ScreenSeqsCommand::screenGroupFile(map<string, string> badSeqNames){
	try {
		ifstream inputGroups;
		m->openInputFile(groupfile, inputGroups);
		string seqName, group;
		map<string, string>::iterator it;
		map<string, string> variables;
		variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(groupfile));
        variables["[extension]"] = m->getExtension(groupfile);
        string goodGroupFile = getOutputFileName("group", variables);
        outputNames.push_back(goodGroupFile);  outputTypes["group"].push_back(goodGroupFile);
		ofstream goodGroupOut;	m->openOutputFile(goodGroupFile, goodGroupOut);
		
		while(!inputGroups.eof()){
			if (m->control_pressed) { goodGroupOut.close(); inputGroups.close(); m->mothurRemove(goodGroupFile); return 0; }

			inputGroups >> seqName; m->gobble(inputGroups); inputGroups >> group; m->gobble(inputGroups);
			it = badSeqNames.find(seqName);
			
			if(it != badSeqNames.end()){
				badSeqNames.erase(it);
			}
			else{
				goodGroupOut << seqName << '\t' << group << endl;
			}
		}
		
		if (m->control_pressed) { goodGroupOut.close();  inputGroups.close(); m->mothurRemove(goodGroupFile);  return 0; }

		//we were unable to remove some of the bad sequences
		if (badSeqNames.size() != 0) {
			for (it = badSeqNames.begin(); it != badSeqNames.end(); it++) {  
				m->mothurOut("Your groupfile does not include the sequence " + it->first + " please correct."); 
				m->mothurOutEndLine();
			}
		}
		
		inputGroups.close();
		goodGroupOut.close();
		
		if (m->control_pressed) { m->mothurRemove(goodGroupFile);   }
		
		return 0;
	
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "screenGroupFile");
		exit(1);
	}
}
//***************************************************************************************************************
int ScreenSeqsCommand::screenCountFile(map<string, string> badSeqNames){
	try {
		ifstream in;
		m->openInputFile(countfile, in);
		map<string, string>::iterator it;
		map<string, string> variables;
		variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(countfile));
        variables["[extension]"] = m->getExtension(countfile);
        string goodCountFile = getOutputFileName("count", variables);
		
        outputNames.push_back(goodCountFile);  outputTypes["count"].push_back(goodCountFile);
		ofstream goodCountOut;	m->openOutputFile(goodCountFile, goodCountOut);
		
        string headers = m->getline(in); m->gobble(in);
        goodCountOut << headers << endl;
        string test = headers; vector<string> pieces = m->splitWhiteSpace(test);
        
        string name, rest; int thisTotal; rest = "";
        while (!in.eof()) {

			if (m->control_pressed) { goodCountOut.close(); in.close(); m->mothurRemove(goodCountFile); return 0; }
            
			in >> name; m->gobble(in); 
            in >> thisTotal; 
            if (pieces.size() > 2) {  rest = m->getline(in); m->gobble(in);  }
            
			it = badSeqNames.find(name);
			
			if(it != badSeqNames.end()){
				badSeqNames.erase(it); 
			}
			else{
				goodCountOut << name << '\t' << thisTotal << '\t' << rest << endl;
			}
		}
		if (m->control_pressed) { goodCountOut.close();  in.close(); m->mothurRemove(goodCountFile);  return 0; }
        
		//we were unable to remove some of the bad sequences
		if (badSeqNames.size() != 0) {
			for (it = badSeqNames.begin(); it != badSeqNames.end(); it++) {  
				m->mothurOut("Your count file does not include the sequence " + it->first + " please correct."); 
				m->mothurOutEndLine();
			}
		}
		
		in.close();
		goodCountOut.close();
        
        //check for groups that have been eliminated
        CountTable ct;
        if (ct.testGroups(goodCountFile)) {
            ct.readTable(goodCountFile, true, false);
            ct.printTable(goodCountFile);
        }
		
		if (m->control_pressed) { m->mothurRemove(goodCountFile);   }
		
		return 0;
        
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "screenCountFile");
		exit(1);
	}
}
//***************************************************************************************************************

int ScreenSeqsCommand::screenTaxonomy(map<string, string> badSeqNames){
	try {
		ifstream input;
		m->openInputFile(taxonomy, input);
		string seqName, tax;
		map<string, string>::iterator it;
        map<string, string> variables;
		variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(taxonomy));
        variables["[extension]"] = m->getExtension(taxonomy);
        string goodTaxFile = getOutputFileName("taxonomy", variables);

		outputNames.push_back(goodTaxFile);  outputTypes["taxonomy"].push_back(goodTaxFile);
		ofstream goodTaxOut;	m->openOutputFile(goodTaxFile, goodTaxOut);
				
		while(!input.eof()){
			if (m->control_pressed) { goodTaxOut.close(); input.close(); m->mothurRemove(goodTaxFile); return 0; }
			
			input >> seqName; m->gobble(input); input >> tax;
			it = badSeqNames.find(seqName);
			
			if(it != badSeqNames.end()){ badSeqNames.erase(it); }
			else{
				goodTaxOut << seqName << '\t' << tax << endl;
			}
			m->gobble(input);
		}
		
		if (m->control_pressed) { goodTaxOut.close(); input.close(); m->mothurRemove(goodTaxFile); return 0; }
		
		//we were unable to remove some of the bad sequences
		if (badSeqNames.size() != 0) {
			for (it = badSeqNames.begin(); it != badSeqNames.end(); it++) {  
				m->mothurOut("Your taxonomy file does not include the sequence " + it->first + " please correct."); 
				m->mothurOutEndLine();
			}
		}
		
		input.close();
		goodTaxOut.close();
		
		if (m->control_pressed) {  m->mothurRemove(goodTaxFile);  return 0; }
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "screenTaxonomy");
		exit(1);
	}
	
}
//***************************************************************************************************************

int ScreenSeqsCommand::screenQual(map<string, string> badSeqNames){
	try {
		ifstream in;
		m->openInputFile(qualfile, in);
		map<string, string>::iterator it;
		map<string, string> variables;
		variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(qualfile));
        variables["[extension]"] = m->getExtension(qualfile);
        string goodQualFile = getOutputFileName("qfile", variables);
		
		outputNames.push_back(goodQualFile);  outputTypes["qfile"].push_back(goodQualFile);
		ofstream goodQual;	m->openOutputFile(goodQualFile, goodQual);
		
		while(!in.eof()){	
			
			if (m->control_pressed) { goodQual.close(); in.close(); m->mothurRemove(goodQualFile); return 0; }

			string saveName = "";
			string name = "";
			string scores = "";
			
			in >> name; 
			
			if (name.length() != 0) { 
				saveName = name.substr(1);
				while (!in.eof())	{	
					char c = in.get(); 
					if (c == 10 || c == 13 || c == -1){	break;	}
					else { name += c; }	
				} 
				m->gobble(in);
			}
			
			while(in){
				char letter= in.get();
				if(letter == '>'){	in.putback(letter);	break;	}
				else{ scores += letter; }
			}
			
			m->gobble(in);
			
			it = badSeqNames.find(saveName);
			
			if(it != badSeqNames.end()){
				badSeqNames.erase(it);
			}else{				
				goodQual << name << endl << scores;
			}
			
			m->gobble(in);
		}
		
		in.close();
		goodQual.close();
		
		//we were unable to remove some of the bad sequences
		if (badSeqNames.size() != 0) {
			for (it = badSeqNames.begin(); it != badSeqNames.end(); it++) {  
				m->mothurOut("Your qual file does not include the sequence " + it->first + " please correct."); 
				m->mothurOutEndLine();
			}
		}
		
		if (m->control_pressed) {  m->mothurRemove(goodQualFile);  return 0; }
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "screenQual");
		exit(1);
	}
	
}
//**********************************************************************************************************************

int ScreenSeqsCommand::driver(linePair filePos, string goodFName, string badAccnosFName, string filename, map<string, string>& badSeqNames){
	try {
		ofstream goodFile;
		m->openOutputFile(goodFName, goodFile);
		
		ofstream badAccnosFile;
		m->openOutputFile(badAccnosFName, badAccnosFile);
		
		ifstream inFASTA;
		m->openInputFile(filename, inFASTA);

		inFASTA.seekg(filePos.start);

		bool done = false;
		int count = 0;
        
		while (!done) {
		
			if (m->control_pressed) {  return 0; }
			
			Sequence currSeq(inFASTA); m->gobble(inFASTA);
			if (currSeq.getName() != "") {
				bool goodSeq = 1;		//	innocent until proven guilty
                string trashCode = "";
                //have the report files found you bad
                map<string, string>::iterator it = badSeqNames.find(currSeq.getName());
                if (it != badSeqNames.end()) { goodSeq = 0;  trashCode = it->second; }  
                
                if (summaryfile == "") { //summaryfile includes these so no need to check again
                    if(startPos != -1 && startPos < currSeq.getStartPos())			{	goodSeq = 0;	trashCode += "start|"; }
                    if(endPos != -1 && endPos > currSeq.getEndPos())				{	goodSeq = 0;	trashCode += "end|";}
                    if(maxAmbig != -1 && maxAmbig <	currSeq.getAmbigBases())		{	goodSeq = 0;	trashCode += "ambig|";}
                    if(maxHomoP != -1 && maxHomoP < currSeq.getLongHomoPolymer())	{	goodSeq = 0;	trashCode += "homop|";}
                    if(minLength != -1 && minLength > currSeq.getNumBases())		{	goodSeq = 0;	trashCode += "<length|";}
                    if(maxLength != -1 && maxLength < currSeq.getNumBases())		{	goodSeq = 0;	trashCode += ">length|";}
                }
                
                if (contigsreport == "") { //contigs report includes this so no need to check again
                    if(maxN != -1 && maxN < currSeq.getNumNs())                     {	goodSeq = 0;	trashCode += "n|"; }
                }
				
				if(goodSeq == 1){
					currSeq.printSequence(goodFile);	
				}else{
					badAccnosFile << currSeq.getName() << '\t' << trashCode.substr(0, trashCode.length()-1) << endl;
					badSeqNames[currSeq.getName()] = trashCode;
				}
                count++;
			}
			
			#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				unsigned long long pos = inFASTA.tellg();
				if ((pos == -1) || (pos >= filePos.end)) { break; }
			#else
				if (inFASTA.eof()) { break; }
			#endif
			
			//report progress
			if((count) % 100 == 0){	m->mothurOutJustToScreen("Processing sequence: " + toString(count)+"\n"); 		}
		}
		//report progress
		if((count) % 100 != 0){	m->mothurOutJustToScreen("Processing sequence: " + toString(count)+"\n"); 	}
		
			
		goodFile.close();
		inFASTA.close();
		badAccnosFile.close();
		
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************
#ifdef USE_MPI
int ScreenSeqsCommand::driverMPI(int start, int num, MPI_File& inMPI, MPI_File& goodFile, MPI_File& badAccnosFile, vector<unsigned long long>& MPIPos, map<string, string>& badSeqNames){
	try {
		string outputString = "";
		MPI_Status statusGood; 
		MPI_Status statusBadAccnos; 
		MPI_Status status; 
		int pid;
		MPI_Comm_rank(MPI_COMM_WORLD, &pid); //find out who we are

		for(int i=0;i<num;i++){
		
			if (m->control_pressed) {  return 0; }
			
			//read next sequence
			int length = MPIPos[start+i+1] - MPIPos[start+i];

			char* buf4 = new char[length];

			MPI_File_read_at(inMPI, MPIPos[start+i], buf4, length, MPI_CHAR, &status);
			
			string tempBuf = buf4;	delete buf4;
			if (tempBuf.length() > length) { tempBuf = tempBuf.substr(0, length);  }
			istringstream iss (tempBuf,istringstream::in);
			
			Sequence currSeq(iss);			
			
			//process seq
			if (currSeq.getName() != "") {
				bool goodSeq = 1;		//	innocent until proven guilty
                string trashCode = "";
                //have the report files found you bad
                map<string, string>::iterator it = badSeqNames.find(currSeq.getName());
                if (it != badSeqNames.end()) { goodSeq = 0;  trashCode = it->second; }  
                
                if (summaryfile == "") { //summaryfile includes these so no need to check again
                    if(startPos != -1 && startPos < currSeq.getStartPos())			{	goodSeq = 0;	trashCode += "start|"; }
                    if(endPos != -1 && endPos > currSeq.getEndPos())				{	goodSeq = 0;	trashCode += "end|";}
                    if(maxAmbig != -1 && maxAmbig <	currSeq.getAmbigBases())		{	goodSeq = 0;	trashCode += "ambig|";}
                    if(maxHomoP != -1 && maxHomoP < currSeq.getLongHomoPolymer())	{	goodSeq = 0;	trashCode += "homop|";}
                    if(minLength != -1 && minLength > currSeq.getNumBases())		{	goodSeq = 0;	trashCode += "<length|";}
                    if(maxLength != -1 && maxLength < currSeq.getNumBases())		{	goodSeq = 0;	trashCode += ">length|";}
                }
                
                if (contigsreport == "") { //contigs report includes this so no need to check again
                    if(maxN != -1 && maxN < currSeq.getNumNs())                     {	goodSeq = 0;	trashCode += "n|"; }
                }
				
                
				if(goodSeq == 1){
					outputString =  ">" + currSeq.getName() + "\n" + currSeq.getAligned() + "\n";
				
					//print good seq
					length = outputString.length();
					char* buf2 = new char[length];
					memcpy(buf2, outputString.c_str(), length);
					
					MPI_File_write_shared(goodFile, buf2, length, MPI_CHAR, &statusGood);
					delete buf2;
				}
				else{

					badSeqNames[currSeq.getName()] = trashCode;
					
					//write to bad accnos file
					outputString = currSeq.getName() + "\t" + trashCode.substr(0, trashCode.length()-1) + "\n";
				
					length = outputString.length();
					char* buf3 = new char[length];
					memcpy(buf3, outputString.c_str(), length);
					
					MPI_File_write_shared(badAccnosFile, buf3, length, MPI_CHAR, &statusBadAccnos);
					delete buf3;
				}
			}
			
			//report progress
			if((i) % 100 == 0){	m->mothurOutJustToScreen("Processing sequence: " + toString(i)+"\n"); 		}
		}
				
		return 1;
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "driverMPI");
		exit(1);
	}
}
#endif
/**************************************************************************************************/

int ScreenSeqsCommand::createProcesses(string goodFileName, string badAccnos, string filename, map<string, string>& badSeqNames) {
	try {
        
        vector<int> processIDS;   
        int process = 1;
		int num = 0;
        bool recalc = false;

#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
				
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driver(lines[process], goodFileName + m->mothurGetpid(process) + ".temp", badAccnos + m->mothurGetpid(process) + ".temp", filename, badSeqNames);
				
				//pass numSeqs to parent
				ofstream out;
				string tempFile = filename + m->mothurGetpid(process) + ".num.temp";
				m->openOutputFile(tempFile, out);
				out << num << endl;
				out.close();
				
				exit(0);
			}else { 
                m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(process) + "\n"); processors = process;
                for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                //wait to die
                for (int i=0;i<processIDS.size();i++) {
                    int temp = processIDS[i];
                    wait(&temp);
                }
                m->control_pressed = false;
                for (int i=0;i<processIDS.size();i++) {
                    m->mothurRemove(filename + (toString(processIDS[i]) + ".num.temp"));
                }
                recalc = true;
                break;
			}
		}
        
        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->control_pressed = false;  for (int i=0;i<processIDS.size();i++) {m->mothurRemove(filename + (toString(processIDS[i]) + ".num.temp"));}processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
            
            //redo file divide
            lines.clear();
            vector<unsigned long long> positions = m->divideFile(filename, processors);
            for (int i = 0; i < (positions.size()-1); i++) {  lines.push_back(linePair(positions[i], positions[(i+1)]));  }
            
            num = 0;
            processIDS.resize(0);
            process = 1;
            
            //loop through and create all the processes you want
            while (process != processors) {
                pid_t pid = fork();
                
                if (pid > 0) {
                    processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
                    process++;
                }else if (pid == 0){
                    num = driver(lines[process], goodFileName + m->mothurGetpid(process) + ".temp", badAccnos + m->mothurGetpid(process) + ".temp", filename, badSeqNames);
                    
                    //pass numSeqs to parent
                    ofstream out;
                    string tempFile = filename + m->mothurGetpid(process) + ".num.temp";
                    m->openOutputFile(tempFile, out);
                    out << num << endl;
                    out.close();
                    
                    exit(0);
                }else { 
                    m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
                    for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                    exit(0);
                }
            }
        }
		
        num = driver(lines[0], goodFileName, badAccnos, filename, badSeqNames);
        
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
		
		for (int i = 0; i < processIDS.size(); i++) {
			ifstream in;
			string tempFile =  filename + toString(processIDS[i]) + ".num.temp";
			m->openInputFile(tempFile, in);
			if (!in.eof()) { int tempNum = 0; in >> tempNum; num += tempNum; }
			in.close(); m->mothurRemove(tempFile);
            
            m->appendFiles((goodFileName + toString(processIDS[i]) + ".temp"), goodFileName);
            m->mothurRemove((goodFileName + toString(processIDS[i]) + ".temp"));
			
            m->appendFiles((badAccnos + toString(processIDS[i]) + ".temp"), badAccnos);
            m->mothurRemove((badAccnos + toString(processIDS[i]) + ".temp"));
		}
		
        //read badSeqs in because root process doesnt know what other "bad" seqs the children found
        ifstream inBad;
        int ableToOpen = m->openInputFile(badAccnos, inBad, "no error");
        
        if (ableToOpen == 0) {
            badSeqNames.clear();
            string tempName, trashCode;
            while (!inBad.eof()) {
                inBad >> tempName >> trashCode; m->gobble(inBad);
                badSeqNames[tempName] = trashCode;
            }
            inBad.close();
        }
#else
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////
		//Windows version shared memory, so be careful when passing variables through the sumScreenData struct. 
		//Above fork() will clone, so memory is separate, but that's not the case with windows, 
		//Taking advantage of shared memory to allow both threads to add info to badSeqNames.
		//////////////////////////////////////////////////////////////////////////////////////////////////////
		
		vector<sumScreenData*> pDataArray; 
		DWORD   dwThreadIdArray[processors-1];
		HANDLE  hThreadArray[processors-1]; 
		
		//Create processor worker threads.
		for( int i=0; i<processors-1; i++ ){
            
            string extension = "";
            if (i!=0) {extension += toString(i) + ".temp"; processIDS.push_back(i); }
            
			// Allocate memory for thread data.
			sumScreenData* tempSum = new sumScreenData(startPos, endPos, maxAmbig, maxHomoP, minLength, maxLength, maxN, badSeqNames, filename, summaryfile, contigsreport, m, lines[i].start, lines[i].end,goodFileName+extension, badAccnos+extension);
			pDataArray.push_back(tempSum);
			
			//default security attributes, thread function name, argument to thread function, use default creation flags, returns the thread identifier
			hThreadArray[i] = CreateThread(NULL, 0, MySumScreenThreadFunction, pDataArray[i], 0, &dwThreadIdArray[i]);   
		}
		
        //do your part
        num = driver(lines[processors-1], (goodFileName+toString(processors-1)+".temp"), (badAccnos+toString(processors-1)+".temp"), filename, badSeqNames);
        processIDS.push_back(processors-1);
        
		//Wait until all threads have terminated.
		WaitForMultipleObjects(processors-1, hThreadArray, TRUE, INFINITE);
		
		//Close all thread handles and free memory allocations.
		for(int i=0; i < pDataArray.size(); i++){
			num += pDataArray[i]->count;
            if (pDataArray[i]->count != pDataArray[i]->end) {
                m->mothurOut("[ERROR]: process " + toString(i) + " only processed " + toString(pDataArray[i]->count) + " of " + toString(pDataArray[i]->end) + " sequences assigned to it, quitting. \n"); m->control_pressed = true; 
            }
            for (map<string, string>::iterator it = pDataArray[i]->badSeqNames.begin(); it != pDataArray[i]->badSeqNames.end(); it++) {	badSeqNames[it->first] = it->second;       }
			CloseHandle(hThreadArray[i]);
			delete pDataArray[i];
		}
        
        for (int i = 0; i < processIDS.size(); i++) {
            m->appendFiles((goodFileName + toString(processIDS[i]) + ".temp"), goodFileName);
            m->mothurRemove((goodFileName + toString(processIDS[i]) + ".temp"));
			
            m->appendFiles((badAccnos + toString(processIDS[i]) + ".temp"), badAccnos);
            m->mothurRemove((badAccnos + toString(processIDS[i]) + ".temp"));
		}

#endif	
        
        return num;
        
	}
	catch(exception& e) {
		m->errorOut(e, "ScreenSeqsCommand", "createProcesses");
		exit(1);
	}
}

//***************************************************************************************************************


