//
//  sffmultiplecommand.cpp
//  Mothur
//
//  Created by Sarah Westcott on 8/14/12.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "sffmultiplecommand.h"



//**********************************************************************************************************************
vector<string> SffMultipleCommand::setParameters(){	
	try {		
		CommandParameter pfile("file", "InputTypes", "", "", "none", "none", "none","fasta-name",false,true,true); parameters.push_back(pfile);
        
        //sffinfo
		CommandParameter ptrim("trim", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(ptrim);
        
        //trim.flows
 		CommandParameter pmaxhomop("maxhomop", "Number", "", "9", "", "", "","",false,false); parameters.push_back(pmaxhomop);
		CommandParameter pmaxflows("maxflows", "Number", "", "450", "", "", "","",false,false); parameters.push_back(pmaxflows);
		CommandParameter pminflows("minflows", "Number", "", "450", "", "", "","",false,false); parameters.push_back(pminflows);
		CommandParameter ppdiffs("pdiffs", "Number", "", "0", "", "", "","",false,false,true); parameters.push_back(ppdiffs);
		CommandParameter pbdiffs("bdiffs", "Number", "", "0", "", "", "","",false,false,true); parameters.push_back(pbdiffs);
        CommandParameter pldiffs("ldiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pldiffs);
		CommandParameter psdiffs("sdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(psdiffs);
        CommandParameter ptdiffs("tdiffs", "Number", "", "0", "", "", "","",false,false); parameters.push_back(ptdiffs);
		CommandParameter psignal("signal", "Number", "", "0.50", "", "", "","",false,false); parameters.push_back(psignal);
		CommandParameter pnoise("noise", "Number", "", "0.70", "", "", "","",false,false); parameters.push_back(pnoise);
		CommandParameter porder("order", "Multiple", "A-B-I", "A", "", "", "","",false,false, true); parameters.push_back(porder);
        //shhh.flows
        CommandParameter plookup("lookup", "InputTypes", "", "", "none", "none", "none","",false,false,true); parameters.push_back(plookup);
		CommandParameter pcutoff("cutoff", "Number", "", "0.01", "", "", "","",false,false); parameters.push_back(pcutoff);
		CommandParameter pmaxiter("maxiter", "Number", "", "1000", "", "", "","",false,false); parameters.push_back(pmaxiter);
        CommandParameter plarge("large", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(plarge);
		CommandParameter psigma("sigma", "Number", "", "60", "", "", "","",false,false); parameters.push_back(psigma);
		CommandParameter pmindelta("mindelta", "Number", "", "0.000001", "", "", "","",false,false); parameters.push_back(pmindelta);
        
        //trim.seqs parameters
        CommandParameter pallfiles("allfiles", "Boolean", "", "t", "", "", "","",false,false); parameters.push_back(pallfiles);
        CommandParameter pflip("flip", "Boolean", "", "F", "", "", "","",false,false,true); parameters.push_back(pflip);
		CommandParameter pmaxambig("maxambig", "Number", "", "-1", "", "", "","",false,false); parameters.push_back(pmaxambig);
		CommandParameter pminlength("minlength", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pminlength);
		CommandParameter pmaxlength("maxlength", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pmaxlength);
		CommandParameter pkeepforward("keepforward", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pkeepforward);
        CommandParameter pkeepfirst("keepfirst", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pkeepfirst);
		CommandParameter premovelast("removelast", "Number", "", "0", "", "", "","",false,false); parameters.push_back(premovelast);

        
        CommandParameter pprocessors("processors", "Number", "", "1", "", "", "","",false,false,true); parameters.push_back(pprocessors);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "SffMultipleCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string SffMultipleCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The sff.multiple command reads a file containing sff filenames and optional oligos filenames. It runs the files through sffinfo, trim.flows, shhh.flows and trim.seqs combining the results.\n";
		helpString += "The sff.multiple command parameters are: ";
        vector<string> parameters = setParameters();
        for (int i = 0; i < parameters.size()-1; i++) {
            helpString += parameters[i] + ", ";
        }
        helpString += parameters[parameters.size()-1] + ".\n";
		helpString += "The file parameter allows you to enter the a file containing the list of sff files and optional oligos files.\n";
        helpString += "The trim parameter allows you to indicate if you would like a sequences and quality scores generated by sffinfo trimmed to the clipQualLeft and clipQualRight values.  Default=True. \n";
        helpString += "The maxambig parameter allows you to set the maximum number of ambiguous bases allowed. The default is -1.\n";
		helpString += "The maxhomop parameter allows you to set a maximum homopolymer length. \n";
		helpString += "The minlength parameter allows you to set and minimum sequence length. \n";
		helpString += "The maxlength parameter allows you to set and maximum sequence length. \n";
		helpString += "The tdiffs parameter is used to specify the total number of differences allowed in the sequence. The default is pdiffs + bdiffs + sdiffs + ldiffs.\n";
		helpString += "The bdiffs parameter is used to specify the number of differences allowed in the barcode. The default is 0.\n";
		helpString += "The pdiffs parameter is used to specify the number of differences allowed in the primer. The default is 0.\n";
        helpString += "The ldiffs parameter is used to specify the number of differences allowed in the linker. The default is 0.\n";
		helpString += "The sdiffs parameter is used to specify the number of differences allowed in the spacer. The default is 0.\n";
		helpString += "The allfiles parameter will create separate group and fasta file for each grouping. The default is F.\n";
		helpString += "The keepforward parameter allows you to indicate whether you want the forward primer removed or not. The default is F, meaning remove the forward primer.\n";
		helpString += "The keepfirst parameter trims the sequence to the first keepfirst number of bases after the barcode or primers are removed, before the sequence is checked to see if it meets the other requirements. \n";
		helpString += "The removelast removes the last removelast number of bases after the barcode or primers are removed, before the sequence is checked to see if it meets the other requirements.\n";
        helpString += "The order parameter options are A, B or I.  Default=A. A = TACG and B = TACGTACGTACGATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGC and I = TACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGC.\n";
		helpString += "Example sff.multiple(file=mySffOligosFile.txt, trim=F).\n";
		helpString += "Note: No spaces between parameter labels (i.e. file), '=' and parameters (i.e.mySffOligosFile.txt).\n";
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "SffMultipleCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string SffMultipleCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "fasta") {  pattern = "[filename],fasta"; } 
        else if (type == "name") {  pattern = "[filename],names"; } 
        else if (type == "group") {  pattern = "[filename],groups"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->control_pressed = true;  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "SffMultipleCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
SffMultipleCommand::SffMultipleCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["fasta"] = tempOutNames;
        outputTypes["name"] = tempOutNames;
        outputTypes["group"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "SffMultipleCommand", "SffMultipleCommand");
		exit(1);
	}
}
//**********************************************************************************************************************

SffMultipleCommand::SffMultipleCommand(string option)  {
	try {
		abort = false; calledHelp = false;  append=false; makeGroup=false;
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			//valid paramters for this command
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
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

			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.validFile(parameters, "outputdir", false);		if (outputDir == "not found"){	outputDir = "";		}
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			inputDir = validParameter.validFile(parameters, "inputdir", false);		
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
                it = parameters.find("file");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["file"] = inputDir + it->second;		}
				}
                
                it = parameters.find("lookup");
				//user has given a template file
				if(it != parameters.end()){ 
					path = m->hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["lookup"] = inputDir + it->second;		}
				}
			}
            
			filename = validParameter.validFile(parameters, "file", true);
            if (filename == "not open") { filename = ""; abort = true; }
            else if (filename == "not found") { filename = "";  }
			
			string temp;
			temp = validParameter.validFile(parameters, "trim", false);					if (temp == "not found"){	temp = "T";				}
			trim = m->isTrue(temp); 
            
            temp = validParameter.validFile(parameters, "minflows", false);	if (temp == "not found") { temp = "450"; }
			m->mothurConvert(temp, minFlows);  
            
			temp = validParameter.validFile(parameters, "maxflows", false);	if (temp == "not found") { temp = "450"; }
			m->mothurConvert(temp, maxFlows);  
            
            temp = validParameter.validFile(parameters, "maxhomop", false);		if (temp == "not found"){	temp = "9";		}
			m->mothurConvert(temp, maxHomoP);  
            
			temp = validParameter.validFile(parameters, "signal", false);		if (temp == "not found"){	temp = "0.50";	}
			m->mothurConvert(temp, signal);  
            
			temp = validParameter.validFile(parameters, "noise", false);		if (temp == "not found"){	temp = "0.70";	}
			m->mothurConvert(temp, noise);  
            
			temp = validParameter.validFile(parameters, "bdiffs", false);		if (temp == "not found"){	temp = "0";		}
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
            
			
			temp = validParameter.validFile(parameters, "processors", false);	if (temp == "not found"){	temp = m->getProcessors();	}
			m->setProcessors(temp);
			m->mothurConvert(temp, processors);
            
			temp = validParameter.validFile(parameters, "order", false);  if (temp == "not found"){ 	temp = "A";	}
            if (temp.length() > 1) {  m->mothurOut("[ERROR]: " + temp + " is not a valid option for order. order options are A, B, or I. A = TACG, B = TACGTACGTACGATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGC, and I = TACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGC.\n");  abort=true;
            }
            else {
                if (toupper(temp[0]) == 'A') {  flowOrder = "A";   }
                else if(toupper(temp[0]) == 'B'){
                    flowOrder = "B";   }
                else if(toupper(temp[0]) == 'I'){
                    flowOrder = "I";   }
                else {
                    m->mothurOut("[ERROR]: " + temp + " is not a valid option for order. order options are A, B, or I. A = TACG, B = TACGTACGTACGATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGC, and I = TACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGC.\n");  abort=true;
                }
            }

            
            temp = validParameter.validFile(parameters, "cutoff", false);	if (temp == "not found"){	temp = "0.01";		}
			m->mothurConvert(temp, cutoff); 
			
			temp = validParameter.validFile(parameters, "mindelta", false);	if (temp == "not found"){	temp = "0.000001";	}
			minDelta = temp; 
            
			temp = validParameter.validFile(parameters, "maxiter", false);	if (temp == "not found"){	temp = "1000";		}
			m->mothurConvert(temp, maxIters); 
            
            temp = validParameter.validFile(parameters, "large", false);	if (temp == "not found"){	temp = "0";		}
			m->mothurConvert(temp, largeSize); 
            if (largeSize != 0) { large = true; }
            else { large = false;  }
            if (largeSize < 0) {  m->mothurOut("The value of the large cannot be negative.\n"); }
            
			temp = validParameter.validFile(parameters, "sigma", false);if (temp == "not found")	{	temp = "60";		}
			m->mothurConvert(temp, sigma); 
            
            temp = validParameter.validFile(parameters, "flip", false);
			if (temp == "not found")    {	flip = 0;	}
			else {  flip = m->isTrue(temp);		}
			
			temp = validParameter.validFile(parameters, "maxambig", false);		if (temp == "not found") { temp = "-1"; }
			m->mothurConvert(temp, maxAmbig);  
                       
			temp = validParameter.validFile(parameters, "minlength", false);	if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, minLength); 
			
			temp = validParameter.validFile(parameters, "maxlength", false);	if (temp == "not found") { temp = "0"; }
			m->mothurConvert(temp, maxLength);
						
			temp = validParameter.validFile(parameters, "keepfirst", false);	if (temp == "not found") { temp = "0"; }
			convert(temp, keepFirst);
            
			temp = validParameter.validFile(parameters, "removelast", false);	if (temp == "not found") { temp = "0"; }
			convert(temp, removeLast);
			
			temp = validParameter.validFile(parameters, "allfiles", false);		if (temp == "not found") { temp = "F"; }
			allFiles = m->isTrue(temp);
            
            temp = validParameter.validFile(parameters, "keepforward", false);		if (temp == "not found") { temp = "F"; }
			keepforward = m->isTrue(temp);
            
            temp = validParameter.validFile(parameters, "lookup", true);
			if (temp == "not found")	{
                string path = m->argv;
                string tempPath = path;
                for (int i = 0; i < path.length(); i++) { tempPath[i] = tolower(path[i]); }
                path = path.substr(0, (tempPath.find_last_of('m')));
                
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
                path += "lookupFiles/";
#else
                path += "lookupFiles\\";
#endif
				lookupFileName = m->getFullPathName(path) + "LookUp_Titanium.pat";
				bool ableToOpen = m->checkLocations(lookupFileName, inputDir);
                if (!ableToOpen) { abort=true; }
			}else if(temp == "not open")	{	
				
				lookupFileName = validParameter.validFile(parameters, "lookup", false);
				
				//if you can't open it its not inputDir, try mothur excutable location
				string exepath = m->argv;
				string tempPath = exepath;
				for (int i = 0; i < exepath.length(); i++) { tempPath[i] = tolower(exepath[i]); }
				exepath = exepath.substr(0, (tempPath.find_last_of('m')));
                
				string tryPath = m->getFullPathName(exepath) + m->getSimpleName(lookupFileName);
				m->mothurOut("Unable to open " + lookupFileName + ". Trying mothur's executable location " + tryPath); m->mothurOutEndLine();
				ifstream in2;
				int ableToOpen = m->openInputFile(tryPath, in2, "noerror");
				in2.close();
				lookupFileName = tryPath;
				
				if (ableToOpen == 1) {  m->mothurOut("Unable to open " + lookupFileName + "."); m->mothurOutEndLine(); abort=true;  }
			}else						{	lookupFileName = temp;	}
		}
	}
	catch(exception& e) {
		m->errorOut(e, "SffMultipleCommand", "SffMultipleCommand");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffMultipleCommand::execute(){
	try {
		if (abort == true) { if (calledHelp) { return 0; }  return 2;	}
		
		vector<string> sffFiles, oligosFiles;
        readFile(sffFiles, oligosFiles);
        
        string thisOutputDir = outputDir;
        if (thisOutputDir == "") { thisOutputDir = m->hasPath(filename); }
        string fileroot = thisOutputDir + m->getRootName(m->getSimpleName(filename));
        map<string, string> variables; 
		variables["[filename]"] = fileroot;
        string fasta = getOutputFileName("fasta",variables);
        string name = getOutputFileName("name",variables);
        string group = getOutputFileName("group",variables);
        
        if (m->control_pressed) { return 0; }
        
        if (sffFiles.size() < processors) { processors = sffFiles.size(); }
        
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#else
        //trim.flows, shhh.flows cannot handle multiple processors for windows.
        processors = 1; m->mothurOut("This command can only use 1 processor on Windows platforms, using 1 processors.\n\n");
#endif
        if (processors == 1) { driver(sffFiles, oligosFiles, 0, sffFiles.size(), fasta, name, group); }
        else { createProcesses(sffFiles, oligosFiles, fasta, name, group); } 
		
		if (m->control_pressed) {  for (int i = 0; i < outputNames.size(); i++) {	m->mothurRemove(outputNames[i]); 	} return 0; }
		
        if (append) { 
            outputNames.push_back(fasta); outputTypes["fasta"].push_back(fasta);
            m->setFastaFile(fasta);
            outputNames.push_back(name); outputTypes["name"].push_back(name);
            m->setNameFile(name);
            if (makeGroup) { outputNames.push_back(group); outputTypes["group"].push_back(group); m->setGroupFile(group); }
        }
        
        m->setProcessors(toString(processors));
        
		//report output filenames
		m->mothurOutEndLine();
		m->mothurOut("Output File Names: "); m->mothurOutEndLine();
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}
		m->mothurOutEndLine();
        
		return 0;
	}
	catch(exception& e) {
		m->errorOut(e, "SffMultipleCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffMultipleCommand::readFile(vector<string>& sffFiles, vector<string>& oligosFiles){
	try {
        
        ifstream in;
        m->openInputFile(filename, in);
        bool allBlank = true;
        bool allFull = true;
        
        string oligos, sff;
        while (!in.eof()) {
            
            if (m->control_pressed) { break; }
            
            in >> sff;
            
            //ignore file pairing
            if(sff[0] == '#'){ while (!in.eof())	{	char c = in.get();  if (c == 10 || c == 13){	break;	}	} m->gobble(in); }
            else { //check for oligos file
                bool ableToOpenSff = m->checkLocations(sff, inputDir);
                
                oligos = "";
            
                // get rest of line in case there is a oligos filename
                while (!in.eof())	{	
                    char c = in.get(); 
                    if (c == 10 || c == 13 || c == -1){	break;	}
                    else if (c == 32 || c == 9){;} //space or tab
                    else { 	oligos += c;  }
                }
                
                if (ableToOpenSff) {
                    sffFiles.push_back(sff);
                    if (oligos != "") {
                        bool ableToOpenOligos = m->checkLocations(oligos, inputDir);
                        if (ableToOpenOligos) {  allBlank = false; }
                        else { m->mothurOut("Can not find " + oligos + ". Ignoring.\n"); oligos = ""; }
                    }
                    if (oligos == "") { allFull = false;  }
                    oligosFiles.push_back(oligos); //will push a blank if there is not an oligos for this sff file
                }else { m->mothurOut("Can not find " + sff + ". Ignoring.\n"); }
            }
            m->gobble(in);
        }
        in.close();
        
        if (allBlank || allFull) { append = true; }
        if (allFull) { makeGroup = true; }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "SffMultipleCommand", "readFile");
		exit(1);
	}
}
//**********************************************************************************************************************
//runs sffinfo, summary.seqs, trim.flows, shhh.flows, trim.seqs, summary.seqs for each sff file.
int SffMultipleCommand::driver(vector<string> sffFiles, vector<string> oligosFiles, int start, int end, string fasta, string name, string group){
    try {
        m->mothurRemove(fasta); m->mothurRemove(name); m->mothurRemove(group);
        int count = 0;
        for (int s = start; s < end; s++) {
            
            string sff = sffFiles[s];
            string oligos = oligosFiles[s];
            
            m->mothurOut("\n>>>>>\tProcessing " + sff + " (file " + toString(s+1) + " of " + toString(sffFiles.size()) + ")\t<<<<<\n");
            
            //run sff.info
            string redirects = "";
            if (inputDir != "")     { redirects += ", inputdir=" + inputDir;    }
            if (outputDir != "")    { redirects += ", outputdir=" + outputDir;  }
            string inputString = "sff=" + sff + ", flow=T";
            if (trim) { inputString += ", trim=T"; }
            if (redirects != "") { inputString += redirects; }
            m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
            m->mothurOut("Running command: sffinfo(" + inputString + ")"); m->mothurOutEndLine(); 
            m->mothurCalling = true;
            
            Command* sffCommand = new SffInfoCommand(inputString);
            sffCommand->execute();
            
            if (m->control_pressed){ break; }
            
            map<string, vector<string> > filenames = sffCommand->getOutputFiles();
            
            delete sffCommand;
            m->mothurCalling = false;
            m->mothurOutEndLine(); 
            
            redirects = "";
            if (outputDir != "")    { redirects += ", outputdir=" + outputDir;  }

            //run summary.seqs on the fasta file
            string fastaFile = "";
            map<string, vector<string> >::iterator it = filenames.find("fasta");
            if (it != filenames.end()) {  if ((it->second).size() != 0) { fastaFile = (it->second)[0];  } }
            else {  m->mothurOut("[ERROR]: sffinfo did not create a fasta file, quitting.\n"); m->control_pressed = true; break;  }
            
            inputString = "fasta=" + fastaFile + ", processors=1";
            if (redirects != "") { inputString += redirects; }
            m->mothurOutEndLine(); 
            m->mothurOut("Running command: summary.seqs(" + inputString + ")"); m->mothurOutEndLine(); 
            m->mothurCalling = true;
            
            Command* summarySeqsCommand = new SeqSummaryCommand(inputString);
            summarySeqsCommand->execute();
            
            if (m->control_pressed){ break; }
            
            map<string, vector<string> > temp = summarySeqsCommand->getOutputFiles();
            mergeOutputFileList(filenames, temp);
            
            delete summarySeqsCommand;
            m->mothurCalling = false;
            
            m->mothurOutEndLine(); 
            
            //run trim.flows on the fasta file
            string flowFile = "";
            it = filenames.find("flow");
            if (it != filenames.end()) {  if ((it->second).size() != 0) { flowFile = (it->second)[0];  } }
            else {  m->mothurOut("[ERROR]: sffinfo did not create a flow file, quitting.\n"); m->control_pressed = true; break;  }
            
            inputString = "flow=" + flowFile;
            if (oligos != "") { inputString += ", oligos=" + oligos; }
            inputString += ", maxhomop=" + toString(maxHomoP) + ", maxflows=" + toString(maxFlows) + ", minflows=" + toString(minFlows);
            inputString += ", pdiffs=" + toString(pdiffs) + ", bdiffs=" + toString(bdiffs) + ", ldiffs=" + toString(ldiffs) + ", sdiffs=" + toString(sdiffs);
            inputString += ", tdiffs=" + toString(tdiffs) + ", signal=" + toString(signal) + ", noise=" + toString(noise) + ", order=" + flowOrder + ", processors=1";
            if (redirects != "") { inputString += redirects; }
            m->mothurOutEndLine(); 
            m->mothurOut("Running command: trim.flows(" + inputString + ")"); m->mothurOutEndLine(); 
            m->mothurCalling = true;
            
            Command* trimFlowCommand = new TrimFlowsCommand(inputString);
            trimFlowCommand->execute();
            
            if (m->control_pressed){ break; }
            
            temp = trimFlowCommand->getOutputFiles();
            mergeOutputFileList(filenames, temp);
            
            delete trimFlowCommand;
            m->mothurCalling = false;
            
            
            string fileFileName = "";
            flowFile = "";
            if (oligos != "") { 
                it = temp.find("file");
                if (it != temp.end()) {  if ((it->second).size() != 0) { fileFileName = (it->second)[0];  } }
                else {  m->mothurOut("[ERROR]: trim.flows did not create a file file, quitting.\n"); m->control_pressed = true; break;  }
            }else {
                vector<string> flowFiles;
                it = temp.find("flow");
                if (it != temp.end()) {  if ((it->second).size() != 0) { flowFiles = (it->second);  } }
                else {  m->mothurOut("[ERROR]: trim.flows did not create a flow file, quitting.\n"); m->control_pressed = true; break;  }
                
                for (int i = 0; i < flowFiles.size(); i++) {
                    string end = flowFiles[i].substr(flowFiles[i].length()-9);
                    if (end == "trim.flow") {
                        flowFile = flowFiles[i]; i+=flowFiles.size(); //if we found the trim.flow file stop looking
                    }
                }
            }
            
            if ((fileFileName == "") && (flowFile == "")) { m->mothurOut("[ERROR]: trim.flows did not create a file file or a trim.flow file, quitting.\n"); m->control_pressed = true; break;  }
            
            if (fileFileName != "") { inputString = "file=" + fileFileName; }
            else { inputString = "flow=" + flowFile; }
            
            inputString += ", lookup=" + lookupFileName + ", cutoff=" + toString(cutoff); + ", maxiters=" + toString(maxIters);
            if (large) { inputString += ", large=" + toString(largeSize); }
            inputString += ", sigma=" +toString(sigma);
            inputString += ", mindelta=" + toString(minDelta);  
            inputString += ", order=" + flowOrder + ", processors=1";
            if (redirects != "") { inputString += redirects; }
            //run shhh.flows
            m->mothurOutEndLine(); 
            m->mothurOut("Running command: shhh.flows(" + inputString + ")"); m->mothurOutEndLine(); 
            m->mothurCalling = true;
            
            Command* shhhFlowCommand = new ShhherCommand(inputString);
            shhhFlowCommand->execute();
            
            if (m->control_pressed){ break; }
            
            temp = shhhFlowCommand->getOutputFiles();
            mergeOutputFileList(filenames, temp);
            
            delete shhhFlowCommand;
            m->mothurCalling = false;
            
            vector<string> fastaFiles;
            vector<string> nameFiles;
            it = temp.find("fasta");
            if (it != temp.end()) {  if ((it->second).size() != 0) { fastaFiles = (it->second);  } }
            else {  m->mothurOut("[ERROR]: shhh.flows did not create a fasta file, quitting.\n"); m->control_pressed = true; break;  }
           
            it = temp.find("name");
            if (it != temp.end()) {  if ((it->second).size() != 0) { nameFiles = (it->second);  } }
            else {  m->mothurOut("[ERROR]: shhh.flows did not create a name file, quitting.\n"); m->control_pressed = true; break;  }
            
            //find fasta and name files with the shortest name.  This is because if there is a composite name it will be the shortest.
            fastaFile = fastaFiles[0];
            for (int i = 1; i < fastaFiles.size(); i++) { if (fastaFiles[i].length() < fastaFile.length()) { fastaFile = fastaFiles[i]; } }
            string nameFile = nameFiles[0];
            for (int i = 1; i < nameFiles.size(); i++) { if (nameFiles[i].length() < nameFile.length()) { nameFile = nameFiles[i]; } }
            
            inputString = "fasta=" + fastaFile + ", name=" + nameFile;
            if (oligos != "") { inputString += ", oligos=" + oligos; }
            if (allFiles) { inputString += ", allfiles=t"; }
            else { inputString += ", allfiles=f";  }
            if (flip) { inputString += ", flip=t"; }
            else { inputString += ", flip=f";  }
            if (keepforward) { inputString += ", keepforward=t"; }
            else { inputString += ", keepforward=f";  }
            
            
            inputString += ", pdiffs=" + toString(pdiffs) + ", bdiffs=" + toString(bdiffs) + ", ldiffs=" + toString(ldiffs) + ", sdiffs=" + toString(sdiffs);
            inputString += ", tdiffs=" + toString(tdiffs) + ", maxambig=" + toString(maxAmbig) + ", minlength=" + toString(minLength) + ", maxlength=" + toString(maxLength);
            if (keepFirst != 0) { inputString += ", keepfirst=" + toString(keepFirst); }
            if (removeLast != 0) { inputString += ", removelast=" + toString(removeLast); }
            inputString += ", processors=1";
            if (redirects != "") { inputString += redirects; }
            //run trim.seqs
            m->mothurOutEndLine(); 
            m->mothurOut("Running command: trim.seqs(" + inputString + ")"); m->mothurOutEndLine(); 
            m->mothurCalling = true;
            
            Command* trimseqsCommand = new TrimSeqsCommand(inputString);
            trimseqsCommand->execute();
            
            if (m->control_pressed){ break; }
            
            temp = trimseqsCommand->getOutputFiles();
            mergeOutputFileList(filenames, temp);
            
            delete trimseqsCommand;
            m->mothurCalling = false;
            
            it = temp.find("fasta");
            if (it != temp.end()) {  if ((it->second).size() != 0) { fastaFiles = (it->second);  } }
            else {  m->mothurOut("[ERROR]: trim.seqs did not create a fasta file, quitting.\n"); m->control_pressed = true; break;  }
            
            for (int i = 0; i < fastaFiles.size(); i++) {
                string end = fastaFiles[i].substr(fastaFiles[i].length()-10);
                if (end == "trim.fasta") {
                    fastaFile = fastaFiles[i]; i+=fastaFiles.size(); //if we found the trim.fasta file stop looking
                }
            }
            
            it = temp.find("name");
            if (it != temp.end()) {  if ((it->second).size() != 0) { nameFiles = (it->second);  } }
            else {  m->mothurOut("[ERROR]: trim.seqs did not create a name file, quitting.\n"); m->control_pressed = true; break;  }
            
            for (int i = 0; i < nameFiles.size(); i++) {
                string end = nameFiles[i].substr(nameFiles[i].length()-10);
                if (end == "trim.names") {
                    nameFile = nameFiles[i]; i+=nameFiles.size(); //if we found the trim.names file stop looking
                }
            }
            
            vector<string> groupFiles;
            string groupFile = "";
            if (makeGroup) {
                it = temp.find("group");
                if (it != temp.end()) {  if ((it->second).size() != 0) { groupFiles = (it->second);  } }
            
                //find group file with the shortest name.  This is because if there is a composite group file it will be the shortest.
                groupFile = groupFiles[0];
                for (int i = 1; i < groupFiles.size(); i++) { if (groupFiles[i].length() < groupFile.length()) { groupFile = groupFiles[i]; } }
            }
            
            inputString = "fasta=" + fastaFile + ", processors=1, name=" + nameFile;
            if (redirects != "") { inputString += redirects; }
            m->mothurOutEndLine(); 
            m->mothurOut("Running command: summary.seqs(" + inputString + ")"); m->mothurOutEndLine(); 
            m->mothurCalling = true;
            
            summarySeqsCommand = new SeqSummaryCommand(inputString);
            summarySeqsCommand->execute();
            
            if (m->control_pressed){ break; }
            
            temp = summarySeqsCommand->getOutputFiles();
            mergeOutputFileList(filenames, temp);
            
            delete summarySeqsCommand;
            m->mothurCalling = false;
            
            m->mothurOutEndLine(); 
            m->mothurOut("/******************************************/"); m->mothurOutEndLine(); 
            
            if (append) {
                m->appendFiles(fastaFile, fasta);
                m->appendFiles(nameFile, name);
                if (makeGroup) { m->appendFiles(groupFile, group);  }
            }
            
            
            for (it = filenames.begin(); it != filenames.end(); it++) {
                for (int i = 0; i < (it->second).size(); i++) {
                    outputNames.push_back((it->second)[i]); outputTypes[it->first].push_back((it->second)[i]);
                }
            }
            count++;
        }
        
        return count;
    }
	catch(exception& e) {
		m->errorOut(e, "SffMultipleCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffMultipleCommand::mergeOutputFileList(map<string, vector<string> >& files, map<string, vector<string> >& temp){
    try {
        map<string, vector<string> >::iterator it;
        for (it = temp.begin(); it != temp.end(); it++) {
            map<string, vector<string> >::iterator it2 = files.find(it->first);
            if (it2 == files.end()) { //we do not already have this type so just add it
                files[it->first] = it->second;
            }else { //merge them
                for (int i = 0; i < (it->second).size(); i++) {
                    files[it->first].push_back((it->second)[i]);
                }
            }
        }
        
        return 0;
    }
	catch(exception& e) {
		m->errorOut(e, "SffMultipleCommand", "mergeOutputFileList");
		exit(1);
	}
}
//**********************************************************************************************************************
int SffMultipleCommand::createProcesses(vector<string> sffFiles, vector<string> oligosFiles, string fasta, string name, string group){
    try {
        vector<int> processIDS;
		int process = 1;
		int num = 0;
        bool recalc = false;
				
		//divide the groups between the processors
		vector<linePair> lines;
        vector<int> numFilesToComplete;
		int numFilesPerProcessor = sffFiles.size() / processors;
		for (int i = 0; i < processors; i++) {
			int startIndex =  i * numFilesPerProcessor;
			int endIndex = (i+1) * numFilesPerProcessor;
			if(i == (processors - 1)){	endIndex = sffFiles.size(); 	}
			lines.push_back(linePair(startIndex, endIndex));
            numFilesToComplete.push_back((endIndex-startIndex));
		}
		
#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)		
		
		//loop through and create all the processes you want
		while (process != processors) {
			pid_t pid = fork();
			
			if (pid > 0) {
				processIDS.push_back(pid);  //create map from line number to pid so you can append files in correct order later
				process++;
			}else if (pid == 0){
				num = driver(sffFiles, oligosFiles, lines[process].start, lines[process].end, fasta + m->mothurGetpid(process) + ".temp", name  + m->mothurGetpid(process) + ".temp", group  + m->mothurGetpid(process) + ".temp");
                
                //pass numSeqs to parent
				ofstream out;
				string tempFile = m->mothurGetpid(process) + ".num.temp";
				m->openOutputFile(tempFile, out);
				out << num << '\t' << outputNames.size() << endl;
                for (int i = 0; i < outputNames.size(); i++) {  out << outputNames[i] << endl;  }
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
                    m->mothurRemove(fasta + (toString(processIDS[i]) + ".temp"));
                    m->mothurRemove(name + (toString(processIDS[i]) + ".temp"));
                    m->mothurRemove(group + (toString(processIDS[i]) + ".temp"));
                }
                recalc = true;
                break;
			}
		}
        
        if (recalc) {
            //test line, also set recalc to true.
            //for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); } for (int i=0;i<processIDS.size();i++) { int temp = processIDS[i]; wait(&temp); } m->control_pressed = false;  for (int i=0;i<processIDS.size();i++) {m->mothurRemove(fasta + (toString(processIDS[i]) + ".temp"));m->mothurRemove(group + (toString(processIDS[i]) + ".temp"));m->mothurRemove(name + (toString(processIDS[i]) + ".temp"));}processors=3; m->mothurOut("[ERROR]: unable to spawn the number of processes you requested, reducing number to " + toString(processors) + "\n");
            
            //redo file divide
            lines.clear();
            numFilesToComplete.clear();
            int numFilesPerProcessor = sffFiles.size() / processors;
            for (int i = 0; i < processors; i++) {
                int startIndex =  i * numFilesPerProcessor;
                int endIndex = (i+1) * numFilesPerProcessor;
                if(i == (processors - 1)){	endIndex = sffFiles.size(); 	}
                lines.push_back(linePair(startIndex, endIndex));
                numFilesToComplete.push_back((endIndex-startIndex));
            }
            
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
                    num = driver(sffFiles, oligosFiles, lines[process].start, lines[process].end, fasta + m->mothurGetpid(process) + ".temp", name  + m->mothurGetpid(process) + ".temp", group  + m->mothurGetpid(process) + ".temp");
                    
                    //pass numSeqs to parent
                    ofstream out;
                    string tempFile = m->mothurGetpid(process) + ".num.temp";
                    m->openOutputFile(tempFile, out);
                    out << num << '\t' << outputNames.size() << endl;
                    for (int i = 0; i < outputNames.size(); i++) {  out << outputNames[i] << endl;  }
                    out.close();
                    
                    exit(0);
                }else { 
                    m->mothurOut("[ERROR]: unable to spawn the necessary processes."); m->mothurOutEndLine(); 
                    for (int i = 0; i < processIDS.size(); i++) { kill (processIDS[i], SIGINT); }
                    exit(0);
                }
            }
        }

		
		//do my part
		num = driver(sffFiles, oligosFiles, lines[0].start, lines[0].end, fasta, name, group);
		
		//force parent to wait until all the processes are done
		for (int i=0;i<processIDS.size();i++) { 
			int temp = processIDS[i];
			wait(&temp);
		}
        
        for (int i=0;i<processIDS.size();i++) { 
            ifstream in;
			string tempFile = toString(processIDS[i]) + ".num.temp";
			m->openInputFile(tempFile, in);
			if (!in.eof()) { 
                int tempNum = 0; int outputNamesSize = 0; 
                in >> tempNum >> outputNamesSize; m->gobble(in);
                for (int j = 0; j < outputNamesSize; j++) {
                    string tempName;
                    in >> tempName; m->gobble(in);
                    outputNames.push_back(tempName);
                }
                if (tempNum != numFilesToComplete[i+1]) {
                    m->mothurOut("[ERROR]: main process expected " + toString(processIDS[i]) + " to complete " + toString(numFilesToComplete[i+1]) + " files, and it only reported completing " + toString(tempNum) + ". This will cause file mismatches.  The flow files may be too large to process with multiple processors. \n");
                }
            }
			in.close(); m->mothurRemove(tempFile);
            
            if (append) {
                m->appendFiles(fasta+toString(processIDS[i])+".temp", fasta);   m->mothurRemove(fasta+toString(processIDS[i])+".temp");
                m->appendFiles(name+toString(processIDS[i])+".temp", name);     m->mothurRemove(name+toString(processIDS[i])+".temp");
                if (makeGroup) { m->appendFiles(group+toString(processIDS[i])+".temp", group);  m->mothurRemove(group+toString(processIDS[i])+".temp"); }
            }
        }
#endif
        return 0;
        
    }
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "createProcesses");
		exit(1);
	}
}
//**********************************************************************************************************************




