/*
 *  trimflowscommand.cpp
 *  Mothur
 *
 *  Created by Pat Schloss on 12/22/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "trimflowscommand.h"
#include "needlemanoverlap.hpp"


//**********************************************************************************************************************
vector<string> TrimFlowsCommand::setParameters(){	
	try {
		CommandParameter pflow("flow", "InputTypes", "", "", "none", "none", "none","flow-file",false,true,true); parameters.push_back(pflow);
		CommandParameter poligos("oligos", "InputTypes", "", "", "none", "none", "none","",false,false,true); parameters.push_back(poligos);
        CommandParameter preorient("checkorient", "Boolean", "", "F", "", "", "","",false,false,true); parameters.push_back(preorient);
		CommandParameter pmaxhomop("maxhomop", "Number", "", "9", "", "", "","",false,false); parameters.push_back(pmaxhomop);
		CommandParameter pmaxflows("maxflows", "Number", "", "450", "", "", "","",false,false); parameters.push_back(pmaxflows);
		CommandParameter pminflows("minflows", "Number", "", "450", "", "", "","",false,false); parameters.push_back(pminflows);
		CommandParameter ppdiffs("pdiffs", "Number", "", "0", "", "", "","",false,false,true); parameters.push_back(ppdiffs);
		CommandParameter pbdiffs("bdiffs", "Number", "", "0", "", "", "","",false,false,true); parameters.push_back(pbdiffs);
        CommandParameter pldiffs("ldiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pldiffs);
		CommandParameter psdiffs("sdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(psdiffs);
        CommandParameter ptdiffs("tdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(ptdiffs);
		CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter psignal("signal", "Number", "", "0.50", "", "", "","",false,false); parameters.push_back(psignal);
		CommandParameter pnoise("noise", "Number", "", "0.70", "", "", "","",false,false); parameters.push_back(pnoise);
		CommandParameter pallfiles("allfiles", "Boolean", "", "t", "", "", "","",false,false); parameters.push_back(pallfiles);
        CommandParameter porder("order", "Multiple", "A-B-I", "A", "", "", "","",false,false, true); parameters.push_back(porder);
		CommandParameter pfasta("fasta", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pfasta);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimFlowsCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string TrimFlowsCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The trim.flows command reads a flowgram file and creates .....\n";
        helpString += "The oligos parameter allows you to provide an oligos file.\n";
        helpString += "The maxhomop parameter allows you to set a maximum homopolymer length. \n";
        helpString += "The tdiffs parameter is used to specify the total number of differences allowed in the sequence. The default is pdiffs + bdiffs + sdiffs + ldiffs.\n";
        helpString += "The checkorient parameter will check look for the reverse compliment of the barcode or primer in the sequence. The default is false.\n";
		helpString += "The bdiffs parameter is used to specify the number of differences allowed in the barcode. The default is 0.\n";
		helpString += "The pdiffs parameter is used to specify the number of differences allowed in the primer. The default is 0.\n";
        helpString += "The ldiffs parameter is used to specify the number of differences allowed in the linker. The default is 0.\n";
		helpString += "The sdiffs parameter is used to specify the number of differences allowed in the spacer. The default is 0.\n";
        helpString += "The order parameter options are A, B or I.  Default=A. A = TACG and B = TACGTACGTACGATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGC and I = TACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGC.\n";
		;
		helpString += "For more details please check out the wiki http://www.mothur.org/wiki/Trim.flows.\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimFlowsCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string TrimFlowsCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "flow") {  pattern = "[filename],[tag],flow"; }
        else if (type == "group") {  pattern = "[filename],flow.groups"; }
        else if (type == "fasta") {  pattern = "[filename],flow.fasta"; } 
        else if (type == "file") {  pattern = "[filename],flow.files"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "TrimFlowsCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************

TrimFlowsCommand::TrimFlowsCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["flow"] = tempOutNames;
        outputTypes["group"] = tempOutNames;
		outputTypes["fasta"] = tempOutNames;
        outputTypes["file"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimFlowsCommand", "TrimFlowsCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

TrimFlowsCommand::TrimFlowsCommand(string option)  {
	try {
		
		abort = false; calledHelp = false;   
		comboStarts = 0;
		
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
				if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["flow"] = tempOutNames;
			outputTypes["fasta"] = tempOutNames;
            outputTypes["file"] = tempOutNames;
            outputTypes["group"] = tempOutNames;
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("flow");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["flow"] = inputDir + it->second;		}
				}
				
				it = parameters.find("oligos");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["oligos"] = inputDir + it->second;		}
				}
				
			}
			
			
			//check for required parameters
			flowFileName = validParameter.validFile(parameters, "flow");
			if (flowFileName == "not found") { 
				flowFileName = current->getFlowFile(); 
				if (flowFileName != "") {  m->mothurOut("Using " + flowFileName + " as input file for the flow parameter."); m->mothurOutEndLine(); }
				else { 
					m->mothurOut("No valid current flow file. You must provide a flow file."); m->mothurOutEndLine(); 
					abort = true;
				} 
			}else if (flowFileName == "not open") { flowFileName = ""; abort = true; }	
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	
				outputDir = "";	
				outputDir += util.hasPath(flowFileName); //if user entered a file with a path then preserve it	
			}
			
			
			//check for optional parameter and set defaults
			// ...at some point should added some additional type checking...
			
			string temp;
			temp = validParameter.valid(parameters, "minflows");	if (temp == "not found") { temp = "450"; }
			util.mothurConvert(temp, minFlows);  

			temp = validParameter.valid(parameters, "maxflows");	if (temp == "not found") { temp = "450"; }
			util.mothurConvert(temp, maxFlows);  
			
			
			temp = validParameter.validFile(parameters, "oligos");
			if (temp == "not found")	{	oligoFileName = "";		}
			else if(temp == "not open")	{	abort = true;			} 
			else						{	oligoFileName = temp;	current->setOligosFile(oligoFileName); }
			
			temp = validParameter.valid(parameters, "fasta");		if (temp == "not found"){	fasta = 0;		}
			else if(util.isTrue(temp))	{	fasta = 1;	}
			
			temp = validParameter.valid(parameters, "maxhomop");		if (temp == "not found"){	temp = "9";		}
			util.mothurConvert(temp, maxHomoP);  

			temp = validParameter.valid(parameters, "signal");		if (temp == "not found"){	temp = "0.50";	}
			util.mothurConvert(temp, signal);  

			temp = validParameter.valid(parameters, "noise");		if (temp == "not found"){	temp = "0.70";	}
			util.mothurConvert(temp, noise);  
	
			temp = validParameter.valid(parameters, "bdiffs");		if (temp == "not found"){	temp = "0";		}
			util.mothurConvert(temp, bdiffs);
			
			temp = validParameter.valid(parameters, "pdiffs");		if (temp == "not found"){	temp = "0";		}
			util.mothurConvert(temp, pdiffs);
			
            temp = validParameter.valid(parameters, "ldiffs");		if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, ldiffs);
            
            temp = validParameter.valid(parameters, "sdiffs");		if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, sdiffs);
			
			temp = validParameter.valid(parameters, "tdiffs");		if (temp == "not found") { int tempTotal = pdiffs + bdiffs + ldiffs + sdiffs;  temp = toString(tempTotal); }
			util.mothurConvert(temp, tdiffs);
			
			if(tdiffs == 0){	tdiffs = bdiffs + pdiffs + ldiffs + sdiffs;	}

			
			temp = validParameter.valid(parameters, "processors");	if (temp == "not found"){	temp = current->getProcessors();	}
			processors = current->setProcessors(temp);
	
			temp = validParameter.valid(parameters, "order");  if (temp == "not found"){ 	temp = "A";	}
            if (temp.length() > 1) {  m->mothurOut("[ERROR]: " + temp + " is not a valid option for order. order options are A, B, or I. A = TACG, B = TACGTACGTACGATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGC, and I = TACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGC.\n");  abort=true;
            }
            else {
                if (toupper(temp[0]) == 'A') {  flowOrder = "TACG";   }
                else if(toupper(temp[0]) == 'B'){
                    flowOrder = "TACGTACGTACGATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGC";   }
                else if(toupper(temp[0]) == 'I'){
                    flowOrder = "TACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGC";   }
                else {
                    m->mothurOut("[ERROR]: " + temp + " is not a valid option for order. order options are A, B, or I. A = TACG, B = TACGTACGTACGATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGC, and I = TACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGC.\n");  abort=true;
                }
            }
            
			if(oligoFileName == "")	{	allFiles = 0;		}
			else					{	allFiles = 1;		}
            
            temp = validParameter.valid(parameters, "checkorient");		if (temp == "not found") { temp = "F"; }
			reorient = util.isTrue(temp);
            
            numBarcodes = 0;
			numFPrimers = 0;
			numRPrimers = 0;
            numLinkers = 0;
            numSpacers = 0;
		}
	}
	catch(exception& e) {
		m->errorOut(e, "TrimFlowsCommand", "TrimFlowsCommand");
		exit(1);
	}
}

//***************************************************************************************************************

int TrimFlowsCommand::execute(){
	try{
		if (abort) { if (calledHelp) { return 0; }  return 2;	}

        map<string, string> variables; 
		variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(flowFileName));
        string fastaFileName = getOutputFileName("fasta",variables);
		if(fasta){ outputNames.push_back(fastaFileName); outputTypes["fasta"].push_back(fastaFileName); }
        
        variables["[tag]"] = "trim";
		string trimFlowFileName = getOutputFileName("flow",variables);
		outputNames.push_back(trimFlowFileName); outputTypes["flow"].push_back(trimFlowFileName);
		
        variables["[tag]"] = "scrap";
		string scrapFlowFileName = getOutputFileName("flow",variables);
		outputNames.push_back(scrapFlowFileName); outputTypes["flow"].push_back(scrapFlowFileName);
        
        createGroup = false;
		if(oligoFileName != ""){   getOligos();	 }
		
        createProcessesCreateTrim(flowFileName, trimFlowFileName, scrapFlowFileName, fastaFileName);
    
		if (m->getControl_pressed()) {  return 0; }			
		
        string flowFilesFileName = getOutputFileName("file",variables);
        outputTypes["file"].push_back(flowFilesFileName);
        outputNames.push_back(flowFilesFileName);
        
		if(allFiles){
            //print group file
            string groupFileName = getOutputFileName("group",variables);
            ofstream out; util.openOutputFile(groupFileName, out);
            for (map<string, string>::iterator it = groupMap.begin(); it != groupMap.end(); it++) {  out << it->first << '\t' << it->second << endl;  } out.close();
            
            //run split.groups command
            string inputString = "flow=" + trimFlowFileName + ", group=" + groupFileName;
            m->mothurOut("/******************************************/\n");
            m->mothurOut("Generating allfiles... Running command: split.groups(" + inputString + ")\n");
            current->setMothurCalling(true);
            
            Command* splitCommand = new SplitGroupCommand(inputString);
            splitCommand->execute();
            
            map<string, vector<string> > filenames = splitCommand->getOutputFiles();
            
            delete splitCommand;
            current->setMothurCalling(false);
            m->mothurOut("/******************************************/\n");

            //print file file
            map<string, vector<string> >::iterator itFiles = filenames.find("flow");
            
            if (itFiles != filenames.end()) {
                ofstream output; util.openOutputFile(flowFilesFileName, output);
                for (int i = 0; i < (itFiles->second).size(); i++) {
                    output << (itFiles->second)[i] << endl;
                    outputNames.push_back((itFiles->second)[i]);
                }
                output.close();
            }else {
                ofstream output; util.openOutputFile(flowFilesFileName, output);
                output << util.getFullPathName(trimFlowFileName) << endl; output.close();
            }
		}else{
            ofstream output; util.openOutputFile(flowFilesFileName, output);
            output << util.getFullPathName(trimFlowFileName) << endl; output.close();
        }
        current->setFileFile(flowFilesFileName);
			
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
		
		return 0;	
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "execute");
		exit(1);
	}
}
/**************************************************************************************************/
struct trimFlowData {
    MothurOut* m;
    string flowFileName, flowOrder;
    OutputWriter* trimFile;
    OutputWriter* scrapFile;
    OutputWriter* fastaFile;
    set<string> badNames;
    unsigned long long lineStart, lineEnd;
    bool pairedOligos, reorient, fasta, createGroup;
    int tdiffs, bdiffs, pdiffs, ldiffs, sdiffs, numFlows, maxHomoP, maxFlows, minFlows;
    float signal, noise;
    long long count;
    vector<string> revPrimer;
    map<string, int> barcodes;
    map<string, int> primers;
    vector<string>  linker;
    vector<string>  spacer;
    vector<string> primerNameVector;
    vector<string> barcodeNameVector;
    map<int, oligosPair> pairedBarcodes;
    map<int, oligosPair> pairedPrimers;
    map<string, string> groupMap;
    Utils util;
    
    trimFlowData(){}
    ~trimFlowData() { }
    trimFlowData(string fn, OutputWriter* tn, OutputWriter* sn, OutputWriter* ffn, bool useFasta, unsigned long long lstart, unsigned long long lend) {
        fasta = useFasta;
        fastaFile = ffn;
        flowFileName = fn;
        trimFile = tn;
        scrapFile = sn;
        lineStart = lstart;
        lineEnd = lend;
        m = MothurOut::getInstance();
    }
    void setOligosOptions(bool cg, int pd, int bd, int ld, int sd, int td, map<string, int> pri, map<string, int> bar, vector<string> revP, vector<string> li, vector<string> spa, map<int, oligosPair> pbr, map<int, oligosPair> ppr, bool po, vector<string> priNameVector, vector<string> barNameVector, bool reo, float sg, float nos, int mhom, string flo, int mxflo, int mnflo, int nmf) {
        createGroup = cg;
        pdiffs = pd;
        bdiffs = bd;
        ldiffs = ld;
        sdiffs = sd;
        tdiffs = td;
        barcodes = bar;
        pairedPrimers = ppr;
        pairedBarcodes = pbr;
        pairedOligos = po;
        primers = pri;
        revPrimer = revP;
        linker = li;
        spacer = spa;
        primerNameVector = priNameVector;
        barcodeNameVector = barNameVector;
        reorient = reo;
        signal = sg;
        noise = nos;
        maxHomoP = mhom;
        flowOrder = flo;
        maxFlows = mxflo;
        minFlows = mnflo;
        numFlows = nmf;
        count = 0;
    }
};
//***************************************************************************************************************

void driverCreateTrim(trimFlowData* params){
	
	try {
        ifstream flowFile; params->util.openInputFile(params->flowFileName, flowFile);
		
		flowFile.seekg(params->lineStart);
		
        if(params->lineStart == 0){ int temp; flowFile >> temp; params->util.gobble(flowFile); }
		
		FlowData flowData(params->numFlows, params->signal, params->noise, params->maxHomoP, params->flowOrder);
		params->count = 0;
		
        int numBarcodes = 0;
        int numLinkers = params->linker.size();
        int numSpacers = params->spacer.size();
        int numFPrimers = 0;
        int numRPrimers = 0;
        TrimOligos* trimOligos = NULL;
        if (params->pairedOligos)   {   trimOligos = new TrimOligos(params->pdiffs, params->bdiffs, 0, 0, params->pairedPrimers, params->pairedBarcodes, false); numBarcodes = params->pairedBarcodes.size(); numFPrimers = params->pairedPrimers.size(); }
        else                {   trimOligos = new TrimOligos(params->pdiffs, params->bdiffs, params->ldiffs, params->sdiffs, params->primers, params->barcodes, params->revPrimer, params->linker, params->spacer); numBarcodes = params->barcodes.size();  numFPrimers = params->primers.size();  numRPrimers = params->revPrimer.size(); }
        
        TrimOligos* rtrimOligos = NULL;
        if (params->reorient) {
            //create reoriented primer and barcode pairs
            map<int, oligosPair> rpairedPrimers, rpairedBarcodes;
            for (map<int, oligosPair>::iterator it = params->pairedPrimers.begin(); it != params->pairedPrimers.end(); it++) {
                oligosPair tempPair(params->util.reverseOligo((it->second).reverse), (params->util.reverseOligo((it->second).forward))); //reversePrimer, rc ForwardPrimer
                rpairedPrimers[it->first] = tempPair;
            }
            for (map<int, oligosPair>::iterator it = params->pairedBarcodes.begin(); it != params->pairedBarcodes.end(); it++) {
                oligosPair tempPair(params->util.reverseOligo((it->second).reverse), (params->util.reverseOligo((it->second).forward))); //reverseBarcode, rc ForwardBarcode
                rpairedBarcodes[it->first] = tempPair;
            }
            int index = rpairedBarcodes.size();
            for (map<string, int>::iterator it = params->barcodes.begin(); it != params->barcodes.end(); it++) {
                oligosPair tempPair("", params->util.reverseOligo((it->first))); //reverseBarcode, rc ForwardBarcode
                rpairedBarcodes[index] = tempPair; index++;
            }
            
            index = rpairedPrimers.size();
            for (map<string, int>::iterator it = params->primers.begin(); it != params->primers.end(); it++) {
                oligosPair tempPair("", params->util.reverseOligo((it->first))); //reverseBarcode, rc ForwardBarcode
                rpairedPrimers[index] = tempPair; index++;
            }
            
            rtrimOligos = new TrimOligos(params->pdiffs, params->bdiffs, 0, 0, rpairedPrimers, rpairedBarcodes, false); numBarcodes = rpairedBarcodes.size();
        }

        bool moreSeqs = 1;
		while(moreSeqs) {
				
			if (params->m->getControl_pressed()) { break; }
			
			int success = 1;
			int currentSeqDiffs = 0;
			string trashCode = "";
            string commentString = "";
			
			flowData.getNext(flowFile); 
			flowData.capFlows(params->maxFlows);
			
			Sequence currSeq = flowData.getSequence();
            //for reorient
            Sequence savedSeq(currSeq.getName(), currSeq.getAligned());
            
			if(!flowData.hasMinFlows(params->minFlows)){	//screen to see if sequence is of a minimum number of flows
				success = 0;
				trashCode += 'l';
			}
            if(!flowData.hasGoodHomoP()){	//screen to see if sequence meets the maximum homopolymer limit
				success = 0;
				trashCode += 'h';
			}

			int primerIndex = 0;
			int barcodeIndex = 0;
			
            if(numLinkers != 0){
                success = trimOligos->stripLinker(currSeq);
                if(success > params->ldiffs)		{	trashCode += 'k';	}
                else{ currentSeqDiffs += success;  }
                
            }
            
            if (params->m->getDebug()) { params->m->mothurOut("[DEBUG]: " + currSeq.getName() + " " + currSeq.getUnaligned() + "\n"); }
            
			if(numBarcodes != 0){
				vector<int> results = trimOligos->stripBarcode(currSeq, barcodeIndex);
                if (params->pairedOligos) {
                    success = results[0] + results[2];
                    commentString += "fbdiffs=" + toString(results[0]) + "(" + trimOligos->getCodeValue(results[1], params->bdiffs) + "), rbdiffs=" + toString(results[2]) + "(" + trimOligos->getCodeValue(results[3], params->bdiffs) + ") ";
                }
                else {
                    success = results[0];
                    commentString += "bdiffs=" + toString(results[0]) + "(" + trimOligos->getCodeValue(results[1], params->bdiffs) + ") ";
                }
				if(success > params->bdiffs)		{	trashCode += 'b';	}
				else{ currentSeqDiffs += success;  }
			}
			
            if(numSpacers != 0){
                success = trimOligos->stripSpacer(currSeq);
                if(success > params->sdiffs)		{	trashCode += 's';	}
                else{ currentSeqDiffs += success;  }
                
            }
            
			if(numFPrimers != 0){
				vector<int> results = trimOligos->stripForward(currSeq, primerIndex);
                if (params->pairedOligos) {
                    success = results[0] + results[2];
                    commentString += "fpdiffs=" + toString(results[0]) + "(" + trimOligos->getCodeValue(results[1], params->pdiffs) + "), rpdiffs=" + toString(results[2]) + "(" + trimOligos->getCodeValue(results[3], params->pdiffs) + ") ";
                }
                else {
                    success = results[0];
                    commentString += "fpdiffs=" + toString(results[0]) + "(" + trimOligos->getCodeValue(results[1], params->pdiffs) + ") ";
                }
				if(success > params->pdiffs)		{	trashCode += 'f';	}
				else{ currentSeqDiffs += success;  }
			}
			
			if(numRPrimers != 0){
                vector<int> results = trimOligos->stripReverse(currSeq);
                success = results[0];
                commentString += "rpdiffs=" + toString(results[0]) + "(" + trimOligos->getCodeValue(results[1], params->pdiffs) + ") ";
                if(success > params->pdiffs)		{	trashCode += 'r';	}
                else{ currentSeqDiffs += success;  }
			}
            
            if (currentSeqDiffs > params->tdiffs)	{	trashCode += 't';   }
            
			if (params->reorient && (trashCode != "")) { //if you failed and want to check the reverse
                int thisSuccess = 0;
                string thisTrashCode = "";
                int thisCurrentSeqsDiffs = 0;
                string thiscommentString = "";
                
                int thisBarcodeIndex = 0;
                int thisPrimerIndex = 0;
               
                if(numBarcodes != 0){
                    vector<int> results = rtrimOligos->stripBarcode(savedSeq, thisBarcodeIndex);
                    if (params->pairedOligos) {
                        thisSuccess = results[0] + results[2];
                        thiscommentString += "fbdiffs=" + toString(results[0]) + "(" + rtrimOligos->getCodeValue(results[1], params->bdiffs) + "), rbdiffs=" + toString(results[2]) + "(" + rtrimOligos->getCodeValue(results[3], params->bdiffs) + ") ";
                    }
                    else {
                        thisSuccess = results[0];
                        thiscommentString += "bdiffs=" + toString(results[0]) + "(" + rtrimOligos->getCodeValue(results[1], params->bdiffs) + ") ";
                    }
                    if(thisSuccess > params->bdiffs)		{ thisTrashCode += "b"; }
                    else{ thisCurrentSeqsDiffs += thisSuccess;  }
                }
               
                if(numFPrimers != 0){
                    vector<int> results = rtrimOligos->stripForward(savedSeq, thisPrimerIndex);
                    if (params->pairedOligos) {
                        thisSuccess = results[0] + results[2];
                        thiscommentString += "fpdiffs=" + toString(results[0]) + "(" + rtrimOligos->getCodeValue(results[1], params->pdiffs) + "), rpdiffs=" + toString(results[2]) + "(" + rtrimOligos->getCodeValue(results[3], params->pdiffs) + ") ";
                    }
                    else {
                        thisSuccess = results[0];
                        thiscommentString += "pdiffs=" + toString(results[0]) + "(" + rtrimOligos->getCodeValue(results[1], params->pdiffs) + ") ";
                    }
                    if(thisSuccess > params->pdiffs)		{ thisTrashCode += "f"; }
                    else{ thisCurrentSeqsDiffs += thisSuccess;  }
                }
                
                if (thisCurrentSeqsDiffs > params->tdiffs)	{	thisTrashCode += 't';   }
                
                if (thisTrashCode == "") {
                    trashCode = thisTrashCode;
                    success = thisSuccess;
                    currentSeqDiffs = thisCurrentSeqsDiffs;
                    commentString = thiscommentString;
                    barcodeIndex = thisBarcodeIndex;
                    primerIndex = thisPrimerIndex;
                    savedSeq.reverseComplement();
                    currSeq.setAligned(savedSeq.getAligned());
                }else { trashCode += "(" + thisTrashCode + ")";  }
            }
            
            currSeq.setComment(commentString);

			if(trashCode.length() == 0){
                
                string thisGroup = "";
                if (params->createGroup) {
                    if(numBarcodes != 0){
                        thisGroup = params->barcodeNameVector[barcodeIndex];
                        if (numFPrimers != 0) {
                            if (params->primerNameVector[primerIndex] != "") {
                                if(thisGroup != "") { thisGroup += "." + params->primerNameVector[primerIndex]; }
                                else                { thisGroup = params->primerNameVector[primerIndex];        }
                            }
                        }
                    }
                }
                
                int pos = thisGroup.find("ignore");
                if (pos == string::npos) {
                    flowData.printFlows(params->trimFile);
                    
                    if(params->fasta)	{ currSeq.printSequence(params->fastaFile);	}
                    
                    if (thisGroup != "") {  params->groupMap[currSeq.getName()] = thisGroup; }
                }
                
			}else{
                params->badNames.insert(currSeq.getName());
                flowData.printFlows(params->scrapFile, trashCode);
            }
            
			params->count++;
            if((params->count) % 10000 == 0){	params->m->mothurOut(toString(params->count)+"\n"); 		}

#if defined NON_WINDOWS
			unsigned long long pos = flowFile.tellg();
			if ((pos == -1) || (pos >= params->lineEnd)) { break; }
#else
			if ((params->count == params->lineEnd) || (flowFile.eof())) { break; }
#endif
		}
        
		//report progress
		if((params->count) % 10000 != 0){	params->m->mothurOut(toString(params->count)+"\n");		}
		
		flowFile.close();
		
        delete trimOligos;
        if (params->reorient) { delete rtrimOligos; }
	}
	catch(exception& e) {
		params->m->errorOut(e, "TrimSeqsCommand", "driverCreateTrim");
		exit(1);
	}
}

//***************************************************************************************************************

int TrimFlowsCommand::getOligos(){
	try {
        bool allBlank = false;
        Oligos oligos; oligos.read(oligoFileName);
        
        if (m->getControl_pressed()) { return 0; } //error in reading oligos
        
        if (oligos.hasPairedBarcodes()) {
            pairedOligos = true;
            pairedPrimers = oligos.getPairedPrimers(); numFPrimers = pairedPrimers.size();
            pairedBarcodes = oligos.getPairedBarcodes(); numBarcodes = pairedBarcodes.size();
        }else {
            pairedOligos = false;
            primers = oligos.getPrimers(); numFPrimers = primers.size();
            barcodes = oligos.getBarcodes(); numBarcodes = barcodes.size();
        }
        
        barcodeNameVector = oligos.getBarcodeNames();
        primerNameVector = oligos.getPrimerNames();
        linker = oligos.getLinkers(); numLinkers = linker.size();
        spacer = oligos.getSpacers(); numSpacers = spacer.size();
        revPrimer = oligos.getReversePrimers(); numRPrimers = revPrimer.size();
        
        vector<string> groupNames = oligos.getGroupNames();
        if (groupNames.size() == 0) { allFiles = 0; allBlank = true;  }
        else { createGroup = true; }
        
        return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimFlowsCommand", "getOligos");
		exit(1);
	}
}
/**************************************************************************************************/
vector<unsigned long long> TrimFlowsCommand::getFlowFileBreaks() {
	try{
		vector<unsigned long long> filePos;
		filePos.push_back(0);
					
		FILE * pFile;
		unsigned long long size;
		
		//get num bytes in file
        flowFileName = util.getFullPathName(flowFileName);
		pFile = fopen (flowFileName.c_str(),"rb");
		if (pFile==NULL) perror ("Error opening file");
		else{
			fseek (pFile, 0, SEEK_END);
			size=ftell (pFile);
			fclose (pFile);
		}
				
		//estimate file breaks
		unsigned long long chunkSize = 0;
		chunkSize = size / processors;

		//file too small to divide by processors
		if (chunkSize == 0)  {  processors = 1;	filePos.push_back(size); return filePos;	}
		
		//for each process seekg to closest file break and search for next '>' char. make that the filebreak
		for (int i = 0; i < processors; i++) {
			unsigned long long spot = (i+1) * chunkSize;
			
			ifstream in;
			util.openInputFile(flowFileName, in);
			in.seekg(spot);
			
			string dummy = util.getline(in);
			
			//there was not another sequence before the end of the file
			unsigned long long sanityPos = in.tellg();
			
//			if (sanityPos == -1) {	break;  }
//			else {  filePos.push_back(newSpot);  }
			if (sanityPos == -1) {	break;  }
			else {  filePos.push_back(sanityPos);  }
			
			in.close();
		}
		
		//save end pos
		filePos.push_back(size);
		
		//sanity check filePos
		for (int i = 0; i < (filePos.size()-1); i++) {
			if (filePos[(i+1)] <= filePos[i]) {  filePos.erase(filePos.begin()+(i+1)); i--; }
		}

		ifstream in;
		util.openInputFile(flowFileName, in);
		in >> numFlows;
		util.gobble(in);
		in.close();
		
		processors = (filePos.size() - 1);
		
		return filePos;	
	}
	catch(exception& e) {
		m->errorOut(e, "TrimSeqsCommand", "getFlowFileBreaks");
		exit(1);
	}
}

/**************************************************************************************************/

int TrimFlowsCommand::createProcessesCreateTrim(string flowFileName, string trimFlowFileName, string scrapFlowFileName, string fastaFileName){

	try {
        time_t start = time(NULL);
        ifstream in; util.openInputFile(flowFileName, in); in >> numFlows; in.close();
        
        vector<linePair> lines;
#if defined NON_WINDOWS
        vector<unsigned long long> flowFilePos = getFlowFileBreaks();
        for (int i = 0; i < (flowFilePos.size()-1); i++) { lines.push_back(linePair(flowFilePos[i], flowFilePos[(i+1)])); }
#else
        
        if (processors == 1) { lines.push_back(linePair(0, -1)); }
        else {
            long long numFlowLines;
            vector<unsigned long long> flowFilePos = util.setFilePosEachLine(flowFileName, numFlowLines);
            
            //figure out how many sequences you have to process
            int numSeqsPerProcessor = numFlowLines / processors;
            
            for (int i = 0; i < processors; i++) {
                int startIndex =  i * numSeqsPerProcessor;
                if(i == (processors - 1)){	numSeqsPerProcessor = numFlowLines - i * numSeqsPerProcessor; 	}
                lines.push_back(linePair(flowFilePos[startIndex], numSeqsPerProcessor));
            }
        }
#endif
        
        //create array of worker threads
        vector<thread*> workerThreads;
        vector<trimFlowData*> data;
        
        ofstream outTrim, outScrap;
        util.openOutputFile(trimFlowFileName, outTrim); outTrim << maxFlows << endl; outTrim.close();
        util.openOutputFile(scrapFlowFileName, outScrap); outScrap << numFlows << endl; outScrap.close();
        
        auto synchronizedOutputTrimFile = std::make_shared<SynchronizedOutputFile>(trimFlowFileName, true); //append
        auto synchronizedOutputScrapFile = std::make_shared<SynchronizedOutputFile>(scrapFlowFileName, true); //append
        auto synchronizedOutputFastaFile = std::make_shared<SynchronizedOutputFile>(fastaFileName);
        
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            OutputWriter* threadTrimWriter = new OutputWriter(synchronizedOutputTrimFile);
            OutputWriter* threadScrapWriter = new OutputWriter(synchronizedOutputScrapFile);
            OutputWriter* threadFastaWriter = NULL;
            
            if (fasta) { threadFastaWriter = new OutputWriter(synchronizedOutputFastaFile); }
            
            trimFlowData* dataBundle = new trimFlowData(flowFileName, threadTrimWriter, threadScrapWriter, threadFastaWriter, fasta, lines[i+1].start, lines[i+1].end);
            dataBundle->setOligosOptions(createGroup, pdiffs, bdiffs, ldiffs, sdiffs, tdiffs, primers, barcodes, revPrimer, linker, spacer, pairedBarcodes, pairedPrimers, pairedOligos,
                                         primerNameVector, barcodeNameVector, reorient, signal, noise, maxHomoP, flowOrder, maxFlows, minFlows, numFlows);
            data.push_back(dataBundle);
            
            workerThreads.push_back(new thread(driverCreateTrim, dataBundle));
        }
        
        OutputWriter* threadTrimWriter = new OutputWriter(synchronizedOutputTrimFile);
        OutputWriter* threadScrapWriter = new OutputWriter(synchronizedOutputScrapFile);
        OutputWriter* threadFastaWriter = NULL;
        
        if (fasta) { threadFastaWriter = new OutputWriter(synchronizedOutputFastaFile); }
        
        trimFlowData* dataBundle = new trimFlowData(flowFileName, threadTrimWriter, threadScrapWriter, threadFastaWriter, fasta, lines[0].start, lines[0].end);
        dataBundle->setOligosOptions(createGroup, pdiffs, bdiffs, ldiffs, sdiffs, tdiffs, primers, barcodes, revPrimer, linker, spacer, pairedBarcodes, pairedPrimers, pairedOligos,
                                     primerNameVector, barcodeNameVector, reorient, signal, noise, maxHomoP, flowOrder, maxFlows, minFlows, numFlows);
        
        driverCreateTrim(dataBundle);
        long long num = dataBundle->count;

        set<string> badNames = dataBundle->badNames;
        groupMap = dataBundle->groupMap;

        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            num += data[i]->count;
            
            delete data[i]->trimFile;
            delete data[i]->scrapFile;
            if (fasta) { delete data[i]->fastaFile; }

            badNames.insert(data[i]->badNames.begin(), data[i]->badNames.end());
            groupMap.insert(data[i]->groupMap.begin(), data[i]->groupMap.end());
            
            delete data[i];
            delete workerThreads[i];
        }
        
       	delete threadTrimWriter;
        delete threadScrapWriter;
        if (fasta) { delete threadFastaWriter; }
        delete dataBundle;
        
        m->mothurOut("It took " + toString(time(NULL) - start) + " secs to trim " + toString(num) + " sequences."); if (m->getDebug()) {   m->mothurOut("Scrapped " + toString(badNames.size()) + ".");  } m->mothurOutEndLine();
        
		return num;
	}
	catch(exception& e) {
		m->errorOut(e, "TrimFlowsCommand", "createProcessesCreateTrim");
		exit(1);
	}
}

//***************************************************************************************************************
