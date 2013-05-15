//
//  makelookupcommand.cpp
//  Mothur
//
//  Created by SarahsWork on 5/14/13.
//  Copyright (c) 2013 Schloss Lab. All rights reserved.
//

#include "makelookupcommand.h"

//**********************************************************************************************************************
vector<string> MakeLookupCommand::setParameters(){
	try {
		CommandParameter ptemplate("reference", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(ptemplate);
        CommandParameter pflow("flow", "InputTypes", "", "", "none", "none", "none","lookup",false,true,true); parameters.push_back(pflow);
        CommandParameter perrors("error", "InputTypes", "", "", "none", "none", "none","none",false,true,true); parameters.push_back(perrors);
        CommandParameter pbarcode("barcode", "String", "", "AACCGTGTC", "", "", "","",false,false); parameters.push_back(pbarcode);
		CommandParameter pkey("key", "String", "", "TCAG", "", "", "","",false,false); parameters.push_back(pkey);
        CommandParameter pthreshold("threshold", "Number", "", "10000", "", "", "","",false,false); parameters.push_back(pthreshold);
        CommandParameter porder("order", "Multiple", "A-B-I", "A", "", "", "","",false,false, true); parameters.push_back(porder);
		CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeLookupCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string MakeLookupCommand::getHelpString(){
	try {
		string helpString = "";
		helpString += "The make.lookup command allows you to create custom lookup files for use with shhh.flows.\n";
		helpString += "The make.lookup command parameters are: reference, flow, error, barcode, key, threshold and order.\n";
		helpString += "The reference file needs to be in the same direction as the flow data and it must start with the forward primer sequence. It is required.\n";
        helpString += "The flow parameter is used to provide the flow data. It is required.\n";
        helpString += "The error parameter is used to provide the error summary. It is required.\n";
        helpString += "The barcode parameter is used to provide the barcode sequence. Default=AACCGTGTC.\n";
        helpString += "The key parameter is used to provide the key sequence. Default=TACG.\n";
        helpString += "The threshold parameter is ....\n";
        helpString += "The order parameter options are A, B or I.  Default=A. A = TACG and B = TACGTACGTACGATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGC and I = TACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGC.\n";
		helpString += "The make.lookup should be in the following format: make.lookup(reference=HMP_MOCK.v53.fasta, flow=H3YD4Z101.mock3.flow_450.flow, error=H3YD4Z101.mock3.flow_450.error.summary, barcode=AACCTGGC)\n";
		helpString += "new(...)\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "MakeLookupCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string MakeLookupCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "lookup") {  pattern = "[filename],lookup"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "MakeLookupCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
MakeLookupCommand::MakeLookupCommand(){
	try {
		abort = true; calledHelp = true;
		setParameters();
        vector<string> tempOutNames;
		outputTypes["lookup"] = tempOutNames; 
    }
	catch(exception& e) {
		m->errorOut(e, "MakeLookupCommand", "MakeLookupCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
MakeLookupCommand::MakeLookupCommand(string option)  {
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
            outputTypes["lookup"] = tempOutNames; 
			
			//if the user changes the input directory command factory will send this info to us in the output parameter
			string inputDir = validParameter.validFile(parameters, "inputdir", false);
			if (inputDir == "not found"){	inputDir = "";		}
			else {
                string path;
				it = parameters.find("flow");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["flow"] = inputDir + it->second;		}
				}
				
				it = parameters.find("error");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["error"] = inputDir + it->second;		}
				}
				
				it = parameters.find("reference");
				//user has given a template file
				if(it != parameters.end()){
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["reference"] = inputDir + it->second;		}
				}
                
            }
                        
			//check for parameters
            errorFileName = validParameter.validFile(parameters, "error", true);
			if (errorFileName == "not open") { errorFileName = ""; abort = true; }
			else if (errorFileName == "not found") { errorFileName = ""; m->mothurOut("[ERROR]: error parameter is required."); m->mothurOutEndLine();  abort = true; }
			
			flowFileName = validParameter.validFile(parameters, "flow", true);
			if (flowFileName == "not open") { flowFileName = ""; abort = true; }
			else if (flowFileName == "not found") { flowFileName = ""; m->mothurOut("[ERROR]: flow parameter is required."); m->mothurOutEndLine();  abort = true; }
			else {   m->setFlowFile(flowFileName);	}
			
			refFastaFileName = validParameter.validFile(parameters, "reference", true);
			if (refFastaFileName == "not open") { abort = true; }
			else if (refFastaFileName == "not found") { refFastaFileName = ""; m->mothurOut("[ERROR]: reference parameter is required."); m->mothurOutEndLine();  abort = true; }
                      
            //if the user changes the output directory command factory will send this info to us in the output parameter
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){
				outputDir = m->hasPath(flowFileName); //if user entered a file with a path then preserve it
			}
            
            string temp = validParameter.validFile(parameters, "threshold", false);	if (temp == "not found"){	temp = "10000";	}
			m->mothurConvert(temp, thresholdCount);
            
            barcodeSequence = validParameter.validFile(parameters, "barcode", false);	if (barcodeSequence == "not found"){	barcodeSequence = "AACCGTGTC";	}
            
            keySequence = validParameter.validFile(parameters, "key", false);	if (keySequence == "not found"){	keySequence = "TCAG";	}
            
            temp = validParameter.validFile(parameters, "order", false);  if (temp == "not found"){ 	temp = "A";	}
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
			
		}
		
	}
	catch(exception& e) {
		m->errorOut(e, "MakeLookupCommand", "MakeLookupCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

int MakeLookupCommand::execute(){
	try {
        
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
        
        cout.setf(ios::fixed, ios::floatfield);
        cout.setf(ios::showpoint);
        
        double gapOpening = 10;
        int maxHomoP = 101;
        vector<vector<double> > penaltyMatrix(maxHomoP);
        for(int i=0;i<maxHomoP;i++){
            penaltyMatrix[i].resize(maxHomoP, 5);
            penaltyMatrix[i][i] = 0;
        }
        
        //Create flows for reference sequences...
        ifstream refFASTA;
        m->openInputFile(refFastaFileName, refFASTA);                        //  *   open reference sequence file
        map<string, vector<double> > refFlowgrams;
        
        while(!refFASTA.eof()){
            if (m->control_pressed) { refFASTA.close(); return 0; }
            Sequence seq(refFASTA);  m->gobble(refFASTA);
            
            string fullSequence = keySequence + barcodeSequence + seq.getAligned(); //  *   concatenate the keySequence, barcodeSequence, and
            //      referenceSequences
            refFlowgrams[seq.getName()] = convertSeqToFlow(fullSequence, flowOrder); //  *   translate concatenated sequences into flowgram
        }
        refFASTA.close();
        
        vector<vector<double> > lookupTable(1000);
        for(int i=0;i<1000;i++){
            lookupTable[i].resize(11, 0);
        }
        
        
        //Loop through each sequence in the flow file and the error summary file.
        ifstream flowFile;
        m->openInputFile(flowFileName, flowFile);
        int numFlows;
        flowFile >> numFlows;
        
        ifstream errorFile;
        m->openInputFile(errorFileName, errorFile);
        m->getline(errorFile); //grab headers
        
        string errorQuery, flowQuery, referenceName, dummy;
        string chimera;
        float intensity;
        
        vector<double> std(11, 0);
        
        while(errorFile && flowFile){
            
            if (m->control_pressed) { errorFile.close(); flowFile.close(); return 0; }
            
            //  * if it's chimeric, chuck it
            errorFile >> errorQuery >> referenceName;
            for(int i=2;i<40;i++){
                errorFile >> dummy;
            }
            errorFile >> chimera;
           
            
            if(chimera == "2"){
                m->getline(flowFile);
            }
            else{
                
                flowFile >> flowQuery >> dummy;
                if(flowQuery != errorQuery){    cout << flowQuery << " != " << errorQuery << endl;  }
                
                vector<double> refFlow = refFlowgrams[referenceName];       //  * compare sequence to its closest reference
                
                vector<double> flowgram(numFlows);
                
                for(int i=0;i<numFlows;i++){
                    flowFile >> intensity;
                    flowgram[i] = intensity;// (int)round(100 * intensity);
                }
                m->gobble(flowFile);
                
                alignFlowGrams(flowgram, refFlow, gapOpening, penaltyMatrix, flowOrder);
                
                if (m->control_pressed) { errorFile.close(); flowFile.close(); return 0; }
                
                for(int i=0;i<flowgram.size();i++){
                    int count = (int)round(100*flowgram[i]);
                    if(count > 1000){count = 999;}
                    if(abs(flowgram[i]-refFlow[i])<=0.50){
                        lookupTable[count][int(refFlow[i])]++;               //  * build table
                        std[int(refFlow[i])] += (100*refFlow[i]-count)*(100*refFlow[i]-count);
                    }
                }
                
            }
            m->gobble(errorFile);
            m->gobble(flowFile);
        }
        errorFile.close(); flowFile.close();
        
        //get probabilities
        vector<int> counts(11, 0);
        int totalCount = 0;
        for(int i=0;i<1000;i++){
            for(int j=0;j<11;j++){
                counts[j] += lookupTable[i][j];
                totalCount += lookupTable[i][j];
            }
        }
        
        int N = 11;
        for(int i=0;i<11;i++){
            if(counts[i] < thresholdCount){ N = i; break; }  //bring back
            std[i] = sqrt(std[i]/(double)(counts[i]));  //bring back
        }
        
        regress(std, N);  //bring back
        
        if (m->control_pressed) { return 0; }
        
        double minProbability = 0.1 / (double)totalCount;
        
        //calculate the negative log probabilities of each intensity given the actual homopolymer length; impute with a guassian when counts are too low
        double sqrtTwoPi = 2.50662827463;//pow(2.0 * 3.14159, 0.5);
        
        for(int i=0;i<1000;i++){
            if (m->control_pressed) { return 0; }
            
            for(int j=0;j<N;j++){
                if(lookupTable[i][j] == 0){
                    lookupTable[i][j] = 1;  //bring back
                }
                lookupTable[i][j] = -log(lookupTable[i][j]/double(counts[j]));  //bring back
            }
            
            for(int j=N;j<11;j++){  //bring back
                double normalProbability = 1.0/((double)std[j] * sqrtTwoPi) * exp(-double(i - j*100)* double(i - j*100)/ double(2*std[j]*std[j]));
                if(normalProbability > minProbability){
                    lookupTable[i][j] = -log(normalProbability);
                }
                else{
                    lookupTable[i][j] = -log(minProbability);
                }
            }
        }
        
        
        //calculate the probability of each homopolymer length
        vector<double> negLogHomoProb(11, 0.00);  //bring back
        for(int i=0;i<N;i++){
            negLogHomoProb[i] = -log(counts[i] / (double)totalCount);
        }
        regress(negLogHomoProb, N);
        
        if (m->control_pressed) { return 0; }
        
        //output data table.  column one is the probability of each homopolymer length
        map<string, string> variables;
        variables["[filename]"] = outputDir + m->getRootName(m->getSimpleName(flowFileName));
		string outputFile = getOutputFileName("lookup",variables);
		outputNames.push_back(outputFile); outputTypes["lookup"].push_back(outputFile);
        
        ofstream lookupFile;
        m->openOutputFile(outputFile, lookupFile);
        lookupFile.precision(8);
        
        for(int j=0;j<11;j++){
            //        lookupFile << counts[j];
            lookupFile << showpoint << negLogHomoProb[j]; //bring back
            for(int i=0;i<1000;i++){
                lookupFile << '\t' << lookupTable[i][j];
            }
            lookupFile << endl;
        }
        lookupFile.close();
        
        m->mothurOut("\nData for homopolymer lengths of " + toString(N) + " and longer were imputed for this analysis\n\n");
         
        if (m->control_pressed) { for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); } return 0; }
		
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "MakeLookupCommand", "execute");
		exit(1);
	}
}
//******************************************************************************************************************************

vector<double> MakeLookupCommand::convertSeqToFlow(string sequence, string order){
    try {
        int seqLength = (int)sequence.length();
        int numFlows = (int)order.length();
        vector<double> flowgram;
        
        int orderIndex = 0;
        int sequenceIndex = 0;
        
        while(orderIndex < numFlows && sequenceIndex < seqLength){
            
            if (m->control_pressed) { return flowgram; }
            
            int homopolymerLength = 1;
            
            char base = sequence[sequenceIndex];
            
            while(base == sequence[sequenceIndex+1] && sequenceIndex < seqLength){
                homopolymerLength++;
                sequenceIndex++;
            }
            
            sequenceIndex++;
            
            for(int i=orderIndex; i<orderIndex+numFlows;i++){
                if(order[i%numFlows] == base){
                    //flowgram[i] = homopolymerLength;
                    orderIndex = i%numFlows;
                    break;
                }else {  flowgram.push_back(0); }
            }
            
            //flowgram[orderIndex] = homopolymerLength;
            flowgram.push_back(homopolymerLength);
            
            orderIndex++;
            orderIndex = orderIndex % numFlows;
        }
        
        return flowgram;
    }
	catch(exception& e) {
		m->errorOut(e, "MakeLookupCommand", "convertSeqToFlow");
		exit(1);
	}
}
//******************************************************************************************************************************

int MakeLookupCommand::alignFlowGrams(vector<double>& flowgram, vector<double>& refFlow, double gapOpening, vector<vector<double> > penaltyMatrix, string flowOrder){
    try {
        int numQueryFlows = (int)flowgram.size();
        int numRefFlows = (int)refFlow.size();
        
            //cout << numQueryFlows << '\t' << numRefFlows << endl;
        
        vector<vector<double> > scoreMatrix(numQueryFlows+1);
        vector<vector<char> > directMatrix(numQueryFlows+1);
        
        for(int i=0;i<=numQueryFlows;i++){
            if (m->control_pressed) { return 0; }
            scoreMatrix[i].resize(numRefFlows+1, 0.00);
            directMatrix[i].resize(numRefFlows+1, 'x');
            
            scoreMatrix[i][0] = i * gapOpening;
            directMatrix[i][0] = 'u';
        }
        
            //cout << numQueryFlows << '\t' << numRefFlows << endl;
        
        
        for(int i=0;i<=numRefFlows;i++){
            scoreMatrix[0][i] = i * gapOpening;
            directMatrix[0][i] = 'l';
        }
        
        for(int i=1;i<=numQueryFlows;i++){
            for(int j=1;j<=numRefFlows;j++){
                if (m->control_pressed) { return 0; }
                double diagonal = 1000000000;
                if(flowOrder[i%flowOrder.length()] == flowOrder[j%flowOrder.length()]){
                    diagonal = scoreMatrix[i-1][j-1] + penaltyMatrix[round(flowgram[i-1])][refFlow[j-1]];
                }
                double up = scoreMatrix[i-1][j] + gapOpening;
                double left = scoreMatrix[i][j-1] + gapOpening;
                
                double minScore = diagonal;
                char direction = 'd';
                
                if(left < diagonal && left < up){
                    minScore = left;
                    direction = 'l';
                }
                else if(up < diagonal && up < left){
                    minScore = up;
                    direction = 'u';
                }
                
                scoreMatrix[i][j] = minScore;
                directMatrix[i][j] = direction;
                
            }
        }
        
        int minRowIndex = numQueryFlows;
        double minRowScore = scoreMatrix[numQueryFlows][numRefFlows];
        for(int i=0;i<numQueryFlows;i++){
            if (m->control_pressed) { return 0; }
            if(scoreMatrix[i][numRefFlows] < minRowScore){
                minRowScore = scoreMatrix[i][numRefFlows];
                minRowIndex = i;
            }
        }
        
        int minColumnIndex = numRefFlows;
        double minColumnScore = scoreMatrix[numQueryFlows][numRefFlows];
        for(int i=0;i<numQueryFlows;i++){
            if (m->control_pressed) { return 0; }
            if(scoreMatrix[numQueryFlows][i] < minColumnScore){
                minColumnScore = scoreMatrix[numQueryFlows][i];
                minColumnIndex = i;
            }
        }
        
        
        int i=minRowIndex;
        int j= minColumnIndex;
        
        vector<double> newFlowgram;
        vector<double> newRefFlowgram;
        
        while(i > 0 && j > 0){
            if (m->control_pressed) { return 0; }
            if(directMatrix[i][j] == 'd'){
                newFlowgram.push_back(flowgram[i-1]);
                newRefFlowgram.push_back(refFlow[j-1]);
                
                i--;
                j--;
            }
            else if(directMatrix[i][j] == 'l'){
                newFlowgram.push_back(0);
                newRefFlowgram.push_back(refFlow[j-1]);
                
                j--;
            }
            else if(directMatrix[i][j] == 'u'){
                newFlowgram.push_back(flowgram[i-1]);
                newRefFlowgram.push_back(0);
                
                i--;
            }
        }
        
        flowgram = newFlowgram;
        refFlow = newRefFlowgram;
        
        return 0;
    
    }
	catch(exception& e) {
		m->errorOut(e, "MakeLookupCommand", "alignFlowGrams");
		exit(1);
	}
}

//******************************************************************************************************************************

int MakeLookupCommand::regress(vector<double>& data, int N){
    try {
        //fit data for larger values of N
        double xMean = 0;
        double yMean = 0;
        
        for(int i=1;i<N;i++){
            if (m->control_pressed) { return 0; }
            xMean += i;
            yMean += data[i];
        }
        xMean /= (N-1);
        yMean /= (N-1);
        
        double numerator = 0;
        double denomenator = 0;
        for(int i=1;i<N;i++){
            if (m->control_pressed) { return 0; }
            numerator += (i-xMean)*(data[i] - yMean);
            denomenator += (i-xMean) * (i-xMean);
        }
        double slope = numerator / denomenator;
        double intercept = yMean - slope * xMean;
        
        for(int i=N;i<11;i++){
            data[i] = intercept + i * slope;
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "MakeLookupCommand", "regress");
		exit(1);
	}
}

//******************************************************************************************************************************

//**********************************************************************************************************************


