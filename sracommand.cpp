//
//  sracommand.cpp
//  Mothur
//
//  Created by SarahsWork on 10/28/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#include "sracommand.h"
#include "sffinfocommand.h"
#include "parsefastaqcommand.h"

//**********************************************************************************************************************
vector<string> SRACommand::setParameters(){
	try {
        CommandParameter psff("sff", "InputTypes", "", "", "sffFastQFile", "sffFastQFile", "none","xml",false,false); parameters.push_back(psff);
        CommandParameter pgroup("group", "InputTypes", "", "", "groupOligos", "none", "none","",false,false); parameters.push_back(pgroup);
        CommandParameter poligos("oligos", "InputTypes", "", "", "groupOligos", "none", "none","",false,false); parameters.push_back(poligos);
        CommandParameter pfile("file", "InputTypes", "", "", "sffFastQFile", "sffFastQFile", "none","xml",false,false); parameters.push_back(pfile);
		CommandParameter pfastq("fastq", "InputTypes", "", "", "sffFastQFile", "sffFastQFile", "none","xml",false,false); parameters.push_back(pfastq);
        //choose only one multiple options
        CommandParameter pplatform("platform", "Multiple", "454-???-???", "454", "", "", "","",false,false); parameters.push_back(pplatform);
        CommandParameter ppdiffs("pdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(ppdiffs);
		CommandParameter pbdiffs("bdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pbdiffs);
        CommandParameter pldiffs("ldiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pldiffs);
		CommandParameter psdiffs("sdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(psdiffs);
        CommandParameter ptdiffs("tdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(ptdiffs);
        
         //every command must have inputdir and outputdir.  This allows mothur users to redirect input and output files.
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SRACommand::getHelpString(){
	try {
		string helpString = "";
		helpString += "The sra command creates the necessary files for a NCBI submission. The xml file and individual sff or fastq files parsed from the original sff or fastq file.\n";
		helpString += "The sra command parameters are: sff, fastq, file, oligos, pdiffs, bdiffs, ldiffs, sdiffs, tdiffs, group.\n";
        helpString += "The sff parameter is used to provide the original sff file.\n";
		helpString += "The fastq parameter is used to provide the original fastq file.\n";
        helpString += "The oligos parameter is used to provide an oligos file to parse your sff or fastq file by.\n";
        helpString += "The group parameter is used to provide the group file to parse your sff or fastq file by.\n";
		helpString += "The file parameter is used to provide a file containing a list of individual fastq or sff files.\n";
        helpString += "The tdiffs parameter is used to specify the total number of differences allowed in the sequence. The default is pdiffs + bdiffs + sdiffs + ldiffs.\n";
		helpString += "The bdiffs parameter is used to specify the number of differences allowed in the barcode. The default is 0.\n";
		helpString += "The pdiffs parameter is used to specify the number of differences allowed in the primer. The default is 0.\n";
        helpString += "The ldiffs parameter is used to specify the number of differences allowed in the linker. The default is 0.\n";
		helpString += "The sdiffs parameter is used to specify the number of differences allowed in the spacer. The default is 0.\n";

		helpString += "The new command should be in the following format: \n";
		helpString += "new(...)\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string SRACommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "xml") {  pattern = "[filename],xml"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "SRACommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
SRACommand::SRACommand(){
	try {
		abort = true; calledHelp = true;
		setParameters();
        vector<string> tempOutNames;
		outputTypes["xml"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "SRACommand");
		exit(1);
	}
}
//**********************************************************************************************************************
SRACommand::SRACommand(string option)  {
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
            outputTypes["xml"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter
			string inputDir = validParameter.validFile(parameters, "inputdir", false);
			if (inputDir == "not found"){	inputDir = "";		}
			else {
            
                string path;
				it = parameters.find("sff");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["sff"] = inputDir + it->second;		}
				}
				
				it = parameters.find("fastq");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["fastq"] = inputDir + it->second;		}
				}
                
                it = parameters.find("file");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["file"] = inputDir + it->second;		}
				}
                
                it = parameters.find("group");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["group"] = inputDir + it->second;		}
				}
                
                it = parameters.find("oligos");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["oligos"] = inputDir + it->second;		}
				}
            }
            
			//check for parameters
            fastqfile = validParameter.validFile(parameters, "fastq", true);
			if (fastqfile == "not open") { fastqfile = "";  abort = true; }
			else if (fastqfile == "not found") { fastqfile = ""; }
			
			sfffile = validParameter.validFile(parameters, "sff", true);
			if (sfffile == "not open") {  sfffile = "";  abort = true; }
			else if (sfffile == "not found") { sfffile = ""; }
            
            file = validParameter.validFile(parameters, "file", true);
			if (file == "not open") {  file = "";  abort = true; }
			else if (file == "not found") { file = ""; }
            
            groupfile = validParameter.validFile(parameters, "group", true);
			if (groupfile == "not open") {  groupfile = "";  abort = true; }
			else if (groupfile == "not found") { groupfile = ""; }
            else {  m->setGroupFile(groupfile); }
            
            oligosfile = validParameter.validFile(parameters, "oligos", true);
			if (oligosfile == "not found")      {	oligosfile = "";	}
			else if(oligosfile == "not open")	{	abort = true;		}
			else {	m->setOligosFile(oligosfile); }
            
            
            file = validParameter.validFile(parameters, "file", true);
			if (file == "not open") {  file = "";  abort = true; }
			else if (file == "not found") { file = ""; }
			
			if ((fastqfile == "") && (sfffile == "") && (sfffile == "")) {
                m->mothurOut("[ERROR]: You must provide a file, sff file or fastq file before you can use the sra command."); m->mothurOutEndLine(); abort = true;
            }
            
            if ((groupfile != "") && (oligosfile != "")) {
                m->mothurOut("[ERROR]: You may not use a group file and an oligos file, only one."); m->mothurOutEndLine(); abort = true;
            }
            
            if ((fastqfile != "") || (sfffile != "")) {
                if ((groupfile == "") && (oligosfile == "")) {
                    oligosfile = m->getOligosFile();
					if (oligosfile != "") {  m->mothurOut("Using " + oligosfile + " as input file for the oligos parameter."); m->mothurOutEndLine(); }
					else {
						groupfile = m->getGroupFile();
                        if (groupfile != "") {  m->mothurOut("Using " + groupfile + " as input file for the group parameter."); m->mothurOutEndLine(); }
                        else {
                            m->mothurOut("[ERROR]: You must provide groupfile or oligos file if splitting a fastq or sff file."); m->mothurOutEndLine(); abort = true;
                        }
					}
                }
            }
			            
            //use only one Mutliple type
			platform = validParameter.validFile(parameters, "platform", false);
			if (platform == "not found") { platform = "454"; }
			
			if ((platform == "454") || (platform == "????") || (platform == "????") || (platform == "????")) { }
			else { m->mothurOut("Not a valid platform option.  Valid platform options are 454, ...."); m->mothurOutEndLine(); abort = true; }
            
            
            string temp = validParameter.validFile(parameters, "bdiffs", false);		if (temp == "not found"){	temp = "0";		}
			m->mothurConvert(temp, bdiffs);
			
			temp = validParameter.validFile(parameters, "pdiffs", false);		if (temp == "not found"){	temp = "0";		}
			m->mothurConvert(temp, pdiffs);
			
            temp = validParameter.validFile(parameters, "ldiffs", false);		if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, ldiffs);
            
            temp = validParameter.validFile(parameters, "sdiffs", false);		if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, sdiffs);
			
			temp = validParameter.validFile(parameters, "tdiffs", false);		if (temp == "not found") { int tempTotal = pdiffs + bdiffs + ldiffs + sdiffs;  temp = toString(tempTotal); }
			m->mothurConvert(temp, tdiffs);
			
			if(tdiffs == 0){	tdiffs = bdiffs + pdiffs + ldiffs + sdiffs;	}
            			
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "SRACommand");
		exit(1);
	}
}
//**********************************************************************************************************************
int SRACommand::execute(){
	try {
		
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        //parse files
        vector<string> filesBySample;
        isSFF = false;
        
        if (file != "")             {       readFile(filesBySample);        }
        else if (sfffile != "")     {       parseSffFile(filesBySample);    }
        else if (fastqfile != "")   {       parseFastqFile(filesBySample);  }
        
        //create xml file
        
		
        //output files created by command
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
        return 0;
		
    }
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "SRACommand");
		exit(1);
	}
}
//**********************************************************************************************************************
int SRACommand::readFile(vector<string>& files){
	try {
        files.clear();
        
        ifstream in;
        m->openInputFile(file, in);
        
        while(!in.eof()) {
            
            if (m->control_pressed) { break; }
            
            string filename;
            in >> filename; m->gobble(in);
            files.push_back(filename);
        }
        in.close();
        
        if (!m->control_pressed) {
            if (files.size() > 0) {
                int pos = files[0].find(".sff");
                if (pos != string::npos) { isSFF = true; } //these files are sff files
            }
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "readFile");
		exit(1);
	}
}
//**********************************************************************************************************************
int SRACommand::parseSffFile(vector<string>& files){
	try {
        isSFF = true;
        //run sffinfo to parse sff file into individual sampled sff files
        string commandString = "sff=" + sfffile;
        if (groupfile != "") { commandString += ", group=" + groupfile; }
        else if (oligosfile != "") {
            commandString += ", oligos=" + oligosfile;
            //add in pdiffs, bdiffs, ldiffs, sdiffs, tdiffs
            if (pdiffs != 0) { commandString += ", pdiffs=" + toString(pdiffs); }
            if (bdiffs != 0) { commandString += ", bdiffs=" + toString(bdiffs); }
            if (ldiffs != 0) { commandString += ", ldiffs=" + toString(ldiffs); }
            if (sdiffs != 0) { commandString += ", sdiffs=" + toString(sdiffs); }
            if (tdiffs != 0) { commandString += ", tdiffs=" + toString(tdiffs); }
        }
        m->mothurOutEndLine();
        m->mothurOut("/******************************************/"); m->mothurOutEndLine();
        m->mothurOut("Running command: sffinfo(" + commandString + ")"); m->mothurOutEndLine();
        m->mothurCalling = true;
        
        Command* sffinfoCommand = new SffInfoCommand(commandString);
        sffinfoCommand->execute();
        
        map<string, vector<string> > filenames = sffinfoCommand->getOutputFiles();
        map<string, vector<string> >::iterator it = filenames.find("sff");
        if (it != filenames.end()) { files = it->second; }
        else { m->control_pressed = true; } // error in sffinfo
        
        delete sffinfoCommand;
        m->mothurCalling = false;
        m->mothurOut("/******************************************/"); m->mothurOutEndLine();
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "readFile");
		exit(1);
	}
}

//**********************************************************************************************************************
int SRACommand::parseFastqFile(vector<string>& files){
	try {
        
        //run sffinfo to parse sff file into individual sampled sff files
        string commandString = "fastq=" + fastqfile;
        if (groupfile != "") { commandString += ", group=" + groupfile; }
        else if (oligosfile != "") {
            commandString += ", oligos=" + oligosfile;
            //add in pdiffs, bdiffs, ldiffs, sdiffs, tdiffs
            if (pdiffs != 0) { commandString += ", pdiffs=" + toString(pdiffs); }
            if (bdiffs != 0) { commandString += ", bdiffs=" + toString(bdiffs); }
            if (ldiffs != 0) { commandString += ", ldiffs=" + toString(ldiffs); }
            if (sdiffs != 0) { commandString += ", sdiffs=" + toString(sdiffs); }
            if (tdiffs != 0) { commandString += ", tdiffs=" + toString(tdiffs); }
        }
        m->mothurOutEndLine();
        m->mothurOut("/******************************************/"); m->mothurOutEndLine();
        m->mothurOut("Running command: fastq.info(" + commandString + ")"); m->mothurOutEndLine();
        m->mothurCalling = true;
        
        Command* fastqinfoCommand = new ParseFastaQCommand(commandString);
        fastqinfoCommand->execute();
        
        map<string, vector<string> > filenames = fastqinfoCommand->getOutputFiles();
        map<string, vector<string> >::iterator it = filenames.find("fastq");
        if (it != filenames.end()) { files = it->second; }
        else { m->control_pressed = true; } // error in sffinfo
        
        delete fastqinfoCommand;
        m->mothurCalling = false;
        m->mothurOut("/******************************************/"); m->mothurOutEndLine();
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "SRACommand", "readFile");
		exit(1);
	}
}
//**********************************************************************************************************************


