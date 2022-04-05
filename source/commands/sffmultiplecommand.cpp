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
        
        vector<string> tempOutNames;
        outputTypes["fasta"] = tempOutNames;
        outputTypes["name"] = tempOutNames;
        outputTypes["group"] = tempOutNames;
        
        abort = false; calledHelp = false;  append=false; makeGroup=false;

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
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "SffMultipleCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
SffMultipleCommand::SffMultipleCommand(string option) : Command()  {
	try {
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string, string> parameters = parser.getParameters();
			
			ValidParameters validParameter;
 			
			
			filename = validParameter.validFile(parameters, "file");
            if (filename == "not open") { filename = ""; abort = true; }
            else if (filename == "not found") { filename = "";  }
			
			string temp;
			temp = validParameter.valid(parameters, "trim");					if (temp == "not found"){	temp = "T";				}
			trim = util.isTrue(temp); 
            
            temp = validParameter.valid(parameters, "minflows");	if (temp == "not found") { temp = "450"; }
			util.mothurConvert(temp, minFlows);  
            
			temp = validParameter.valid(parameters, "maxflows");	if (temp == "not found") { temp = "450"; }
			util.mothurConvert(temp, maxFlows);  
            
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
                if (toupper(temp[0]) == 'A') {  flowOrder = "A";   }
                else if(toupper(temp[0]) == 'B'){
                    flowOrder = "B";   }
                else if(toupper(temp[0]) == 'I'){
                    flowOrder = "I";   }
                else {
                    m->mothurOut("[ERROR]: " + temp + " is not a valid option for order. order options are A, B, or I. A = TACG, B = TACGTACGTACGATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATAGATCGCATGACGATCGCATATCGTCAGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGTAGTCGAGCATCATCTGACGCAGTACGTGCATGATCTCAGTCAGCAGCTATGTCAGTGCATGCATAGATCGCATGACGATCGCATATCGTCAGTGCAGTGACTGATCGTCATCAGCTAGCATCGACTGCATGATCTCAGTCAGCAGC, and I = TACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGC.\n");  abort=true;
                }
            }

            
            temp = validParameter.valid(parameters, "cutoff");	if (temp == "not found"){	temp = "0.01";		}
			util.mothurConvert(temp, cutoff); 
			
			temp = validParameter.valid(parameters, "mindelta");	if (temp == "not found"){	temp = "0.000001";	}
			minDelta = temp; 
            
			temp = validParameter.valid(parameters, "maxiter");	if (temp == "not found"){	temp = "1000";		}
			util.mothurConvert(temp, maxIters); 
            
            temp = validParameter.valid(parameters, "large");	if (temp == "not found"){	temp = "0";		}
			util.mothurConvert(temp, largeSize); 
            if (largeSize != 0) { large = true; }
            else { large = false;  }
            if (largeSize < 0) {  m->mothurOut("The value of the large cannot be negative.\n"); }
            
			temp = validParameter.valid(parameters, "sigma");if (temp == "not found")	{	temp = "60";		}
			util.mothurConvert(temp, sigma); 
            
            temp = validParameter.valid(parameters, "flip");
			if (temp == "not found")    {	flip = 0;	}
			else {  flip = util.isTrue(temp);		}
			
			temp = validParameter.valid(parameters, "maxambig");		if (temp == "not found") { temp = "-1"; }
			util.mothurConvert(temp, maxAmbig);  
                       
			temp = validParameter.valid(parameters, "minlength");	if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, minLength); 
			
			temp = validParameter.valid(parameters, "maxlength");	if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, maxLength);
						
			temp = validParameter.valid(parameters, "keepfirst");	if (temp == "not found") { temp = "0"; }
			convert(temp, keepFirst);
            
			temp = validParameter.valid(parameters, "removelast");	if (temp == "not found") { temp = "0"; }
			convert(temp, removeLast);
			
			temp = validParameter.valid(parameters, "allfiles");		if (temp == "not found") { temp = "F"; }
			allFiles = util.isTrue(temp);
            
            temp = validParameter.valid(parameters, "keepforward");		if (temp == "not found") { temp = "F"; }
			keepforward = util.isTrue(temp);
            
            temp = validParameter.validFile(parameters, "lookup");
			if (temp == "not found")	{
                string path = current->getProgramPath();
                //string tempPath = path;
                //for (int i = 0; i < path.length(); i++) { tempPath[i] = tolower(path[i]); }
                //path = path.substr(0, (tempPath.find_last_of('m')));
                
#if defined NON_WINDOWS
                path += "lookupFiles/";
#else
                path += "lookupFiles\\";
#endif
				lookupFileName = util.getFullPathName(path) + "LookUp_Titanium.pat";
				bool ableToOpen = util.checkLocations(lookupFileName, current->getLocations());
                if (!ableToOpen) { abort=true; }
			}else if(temp == "not open")	{	
				
				lookupFileName = validParameter.validPath(parameters, "lookup");
				
				//if you can't open it its not inputDir, try mothur excutable location
				string exepath = current->getProgramPath();
				//string tempPath = exepath;
				//for (int i = 0; i < exepath.length(); i++) { tempPath[i] = tolower(exepath[i]); }
				//exepath = exepath.substr(0, (tempPath.find_last_of('m')));
                
				string tryPath = util.getFullPathName(exepath) + util.getSimpleName(lookupFileName);
				m->mothurOut("Unable to open " + lookupFileName + ". Trying mothur's executable location " + tryPath); m->mothurOutEndLine();
				ifstream in2;
				bool ableToOpen = util.openInputFile(tryPath, in2, "noerror");
				in2.close();
				lookupFileName = tryPath;
				
				if (!ableToOpen) {  m->mothurOut("Unable to open " + lookupFileName + ".\n");  abort=true;  }
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
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		vector<string> sffFiles, oligosFiles;
        readFile(sffFiles, oligosFiles);
        
        string thisOutputDir = outputdir;
        if (thisOutputDir == "") { thisOutputDir = util.hasPath(filename); }
        string fileroot = thisOutputDir + util.getRootName(util.getSimpleName(filename));
        map<string, string> variables; 
		variables["[filename]"] = fileroot;
        string fasta = getOutputFileName("fasta",variables);
        string name = getOutputFileName("name",variables);
        string group = getOutputFileName("group",variables);
        
        if (m->getControl_pressed()) { return 0; }
        
        if (sffFiles.size() < processors) { processors = sffFiles.size(); m->mothurOut("Reducing processors to " + toString(sffFiles.size()) + ".\n"); }
    
        createProcesses(sffFiles, oligosFiles, fasta, name, group);
		
		if (m->getControl_pressed()) {  for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]); 	} return 0; }
		
        if (append) { 
            outputNames.push_back(fasta); outputTypes["fasta"].push_back(fasta);
            current->setFastaFile(fasta);
            outputNames.push_back(name); outputTypes["name"].push_back(name);
            current->setNameFile(name);
            if (makeGroup) { outputNames.push_back(group); outputTypes["group"].push_back(group); current->setGroupFile(group); }
        }
        
        current->setProcessors(toString(processors));
        
		//report output filenames
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
        
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
        
        ifstream in; util.openInputFile(filename, in);
        bool allBlank = true;  bool allFull = true;
        
        string oligos, sff;
        while (!in.eof()) {
            
            if (m->getControl_pressed()) { break; }
            
            in >> sff;
            
            //ignore file pairing
            if(sff[0] == '#'){ while (!in.eof())	{	char c = in.get();  if (c == 10 || c == 13){	break;	}	} gobble(in); }
            else { //check for oligos file
                bool ableToOpenSff = util.checkLocations(sff, current->getLocations());
                
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
                        bool ableToOpenOligos = util.checkLocations(oligos, current->getLocations());
                        if (ableToOpenOligos) {  allBlank = false; }
                        else { m->mothurOut("Can not find " + oligos + ". Ignoring.\n"); oligos = ""; }
                    }
                    if (oligos == "") { allFull = false;  }
                    oligosFiles.push_back(oligos); //will push a blank if there is not an oligos for this sff file
                }else { m->mothurOut("Can not find " + sff + ". Ignoring.\n"); }
            }
            gobble(in);
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
void mergeOutputFileList(map<string, vector<string> >& files, map<string, vector<string> >& temp){
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
}
/**************************************************************************************************/
struct sffMultipleData {
    string fasta, name, group;
    vector<string> sffFiles, oligosFiles;
    int start, end;
    MothurOut* m;
    Utils util;
    int count;
    
    string flowOrder, lookupFileName, minDelta;
    bool trim, large, flip, allFiles, keepforward, append, makeGroup;
    int maxFlows, minFlows, minLength, maxLength, maxHomoP, tdiffs, bdiffs, pdiffs, sdiffs, ldiffs;
    int maxIters, largeSize;
    float signal, noise, cutoff, sigma;
    int keepFirst, removeLast, maxAmbig;
    
    vector<string> outputNames;
    map<string, vector<string> > outputTypes;
    
    sffMultipleData(){}
    sffMultipleData(vector<string> sFiles, vector<string> oFiles, string fa, string nm, string grp, int st, int en) {
        sffFiles = sFiles;
        oligosFiles = oFiles;
        fasta = fa;
        name = nm;
        group = grp;
        start = st;
        end = en;
        m = MothurOut::getInstance();
        count = 0;
    }
    
    void setVariables(bool tr, bool lg, bool alf, bool flp, bool kpfo, bool mkg, bool app, int lgs, int bd, int td, int pd, int sd, int ld, int mxf, int mnf, int mnl, int mxl, int mxh, int mxi, int kpf, int rml, float sgn, float n, float cu, float sig, string fo, string lkf, string mnd) {
        trim = tr;
        large = lg;
        allFiles = alf;
        keepforward = kpfo;
        flip = flp;
        makeGroup = mkg;
        append = app;
        largeSize = lgs;
        tdiffs = td;
        bdiffs = bd;
        pdiffs = pd;
        sdiffs = sd;
        ldiffs = ld;
        maxFlows = mxf;
        minFlows = mnf;
        minLength = mnl;
        maxLength = mxl;
        maxHomoP = mxh;
        maxIters = mxi;
        signal = sgn;
        noise = n;
        cutoff = cu;
        sigma = sig;
        flowOrder = fo;
        lookupFileName = lkf;
        minDelta = mnd;
        keepFirst = kpf;
        removeLast = rml;
    }
};
//**********************************************************************************************************************
//runs sffinfo, summary.seqs, trim.flows, shhh.flows, trim.seqs, summary.seqs for each sff file.
void driverSFFMultiple(sffMultipleData* params){
    try {
        params->util.mothurRemove(params->fasta); params->util.mothurRemove(params->name); params->util.mothurRemove(params->group);
        params->count = 0;
        for (int s = params->start; s < params->end; s++) {
            
            string sff = params->sffFiles[s];
            string oligos = params->oligosFiles[s];
            
            params->m->mothurOut("\n>>>>>\tProcessing " + sff + " (file " + toString(s+1) + " of " + toString(params->sffFiles.size()) + ")\t<<<<<\n");
            
            //run sff.info
            string inputString = "sff=" + sff + ", flow=T";
            if (params->trim) { inputString += ", trim=T"; }
        
            params->m->mothurOut("/******************************************/\n");
            params->m->mothurOut("Running command: sffinfo(" + inputString + ")\n");
            
            Command* sffCommand = new SffInfoCommand(inputString);
            sffCommand->execute();
            
            if (params->m->getControl_pressed()){ break; }
            
            map<string, vector<string> > filenames = sffCommand->getOutputFiles();
            delete sffCommand;
            params->m->mothurOutEndLine();

            //run summary.seqs on the fasta file
            string fastaFile = "";
            map<string, vector<string> >::iterator it = filenames.find("fasta");
            if (it != filenames.end()) {  if ((it->second).size() != 0) { fastaFile = (it->second)[0];  } }
            else {  params->m->mothurOut("[ERROR]: sffinfo did not create a fasta file, quitting.\n"); params->m->setControl_pressed(true); break;  }
            
            inputString = "fasta=" + fastaFile + ", processors=1";
            params->m->mothurOut("\nRunning command: summary.seqs(" + inputString + ")\n");
            
            Command* summarySeqsCommand = new SeqSummaryCommand(inputString);
            summarySeqsCommand->execute();
            
            if (params->m->getControl_pressed()){ break; }
            
            map<string, vector<string> > temp = summarySeqsCommand->getOutputFiles();
            mergeOutputFileList(filenames, temp);
            
            delete summarySeqsCommand;
            params->m->mothurOutEndLine();
            
            //run trim.flows on the fasta file
            string flowFile = "";
            it = filenames.find("flow");
            if (it != filenames.end()) {  if ((it->second).size() != 0) { flowFile = (it->second)[0];  } }
            else {  params->m->mothurOut("[ERROR]: sffinfo did not create a flow file, quitting.\n"); params->m->setControl_pressed(true); break;  }
            
            inputString = "flow=" + flowFile;
            if (oligos != "") { inputString += ", oligos=" + oligos; }
            inputString += ", maxhomop=" + toString(params->maxHomoP) + ", maxflows=" + toString(params->maxFlows) + ", minflows=" + toString(params->minFlows);
            inputString += ", pdiffs=" + toString(params->pdiffs) + ", bdiffs=" + toString(params->bdiffs) + ", ldiffs=" + toString(params->ldiffs) + ", sdiffs=" + toString(params->sdiffs);
            inputString += ", tdiffs=" + toString(params->tdiffs) + ", signal=" + toString(params->signal) + ", noise=" + toString(params->noise) + ", order=" + params->flowOrder + ", processors=1";
            params->m->mothurOut("\nRunning command: trim.flows(" + inputString + ")\n");
           
            Command* trimFlowCommand = new TrimFlowsCommand(inputString);
            trimFlowCommand->execute();
            
            if (params->m->getControl_pressed()){ break; }
            
            temp = trimFlowCommand->getOutputFiles();
            mergeOutputFileList(filenames, temp);
            
            delete trimFlowCommand;

            string fileFileName = "";
            flowFile = "";
            if (oligos != "") { 
                it = temp.find("file");
                if (it != temp.end()) {  if ((it->second).size() != 0) { fileFileName = (it->second)[0];  } }
                else {  params->m->mothurOut("[ERROR]: trim.flows did not create a file file, quitting.\n"); params->m->setControl_pressed(true); break;  }
            }else {
                vector<string> flowFiles;
                it = temp.find("flow");
                if (it != temp.end()) {  if ((it->second).size() != 0) { flowFiles = (it->second);  } }
                else {  params->m->mothurOut("[ERROR]: trim.flows did not create a flow file, quitting.\n"); params->m->setControl_pressed(true); break;  }
                
                for (int i = 0; i < flowFiles.size(); i++) {
                    string end = flowFiles[i].substr(flowFiles[i].length()-9);
                    if (end == "trim.flow") {
                        flowFile = flowFiles[i]; i+=flowFiles.size(); //if we found the trim.flow file stop looking
                    }
                }
            }
            
            if ((fileFileName == "") && (flowFile == "")) { params->m->mothurOut("[ERROR]: trim.flows did not create a file file or a trim.flow file, quitting.\n"); params->m->setControl_pressed(true); break;  }
            
            if (fileFileName != "") { inputString = "file=" + fileFileName; }
            else { inputString = "flow=" + flowFile; }
            
            inputString += ", lookup=" + params->lookupFileName + ", cutoff=" + toString(params->cutoff); + ", maxiters=" + toString(params->maxIters);
            if (params->large) { inputString += ", large=" + toString(params->largeSize); }
            inputString += ", sigma=" +toString(params->sigma);
            inputString += ", mindelta=" + toString(params->minDelta);
            inputString += ", order=" + params->flowOrder;
            //run shhh.flows
            params->m->mothurOut("\nRunning command: shhh.flows(" + inputString + ")\n");
            
            Command* shhhFlowCommand = new ShhherCommand(inputString);
            shhhFlowCommand->execute();
            
            if (params->m->getControl_pressed()){ break; }
            
            temp = shhhFlowCommand->getOutputFiles();
            mergeOutputFileList(filenames, temp);
            
            delete shhhFlowCommand;
            
            vector<string> fastaFiles;
            vector<string> nameFiles;
            it = temp.find("fasta");
            if (it != temp.end()) {  if ((it->second).size() != 0) { fastaFiles = (it->second);  } }
            else {  params->m->mothurOut("[ERROR]: shhh.flows did not create a fasta file, quitting.\n"); params->m->setControl_pressed(true); break;  }
           
            it = temp.find("name");
            if (it != temp.end()) {  if ((it->second).size() != 0) { nameFiles = (it->second);  } }
            else {  params->m->mothurOut("[ERROR]: shhh.flows did not create a name file, quitting.\n"); params->m->setControl_pressed(true); break;  }
            
            //find fasta and name files with the shortest name.  This is because if there is a composite name it will be the shortest.
            fastaFile = fastaFiles[0];
            for (int i = 1; i < fastaFiles.size(); i++) { if (fastaFiles[i].length() < fastaFile.length()) { fastaFile = fastaFiles[i]; } }
            string nameFile = nameFiles[0];
            for (int i = 1; i < nameFiles.size(); i++) { if (nameFiles[i].length() < nameFile.length()) { nameFile = nameFiles[i]; } }
            
            inputString = "fasta=" + fastaFile + ", name=" + nameFile;
            if (oligos != "") { inputString += ", oligos=" + oligos; }
            if (params->allFiles) { inputString += ", allfiles=t"; }
            else { inputString += ", allfiles=f";  }
            if (params->flip) { inputString += ", flip=t"; }
            else { inputString += ", flip=f";  }
            if (params->keepforward) { inputString += ", keepforward=t"; }
            else { inputString += ", keepforward=f";  }
            
            
            inputString += ", pdiffs=" + toString(params->pdiffs) + ", bdiffs=" + toString(params->bdiffs) + ", ldiffs=" + toString(params->ldiffs) + ", sdiffs=" + toString(params->sdiffs);
            inputString += ", tdiffs=" + toString(params->tdiffs) + ", maxambig=" + toString(params->maxAmbig) + ", minlength=" + toString(params->minLength) + ", maxlength=" + toString(params->maxLength);
            if (params->keepFirst != 0) { inputString += ", keepfirst=" + toString(params->keepFirst); }
            if (params->removeLast != 0) { inputString += ", removelast=" + toString(params->removeLast); }
            inputString += ", processors=1";
            //run trim.seqs
            params->m->mothurOut("\nRunning command: trim.seqs(" + inputString + ")\n");
            
            Command* trimseqsCommand = new TrimSeqsCommand(inputString);
            trimseqsCommand->execute();
            
            if (params->m->getControl_pressed()){ break; }
            
            temp = trimseqsCommand->getOutputFiles();
            mergeOutputFileList(filenames, temp);
            
            delete trimseqsCommand;
            
            it = temp.find("fasta");
            if (it != temp.end()) {  if ((it->second).size() != 0) { fastaFiles = (it->second);  } }
            else {  params->m->mothurOut("[ERROR]: trim.seqs did not create a fasta file, quitting.\n"); params->m->setControl_pressed(true); break;  }
            
            for (int i = 0; i < fastaFiles.size(); i++) {
                string end = fastaFiles[i].substr(fastaFiles[i].length()-10);
                if (end == "trim.fasta") {
                    fastaFile = fastaFiles[i]; i+=fastaFiles.size(); //if we found the trim.fasta file stop looking
                }
            }
            
            it = temp.find("name");
            if (it != temp.end()) {  if ((it->second).size() != 0) { nameFiles = (it->second);  } }
            else {  params->m->mothurOut("[ERROR]: trim.seqs did not create a name file, quitting.\n"); params->m->setControl_pressed(true); break;  }
            
            for (int i = 0; i < nameFiles.size(); i++) {
                string end = nameFiles[i].substr(nameFiles[i].length()-10);
                if (end == "trim.names") {
                    nameFile = nameFiles[i]; i+=nameFiles.size(); //if we found the trim.names file stop looking
                }
            }
            
            vector<string> groupFiles;
            string groupFile = "";
            if (params->makeGroup) {
                it = temp.find("group");
                if (it != temp.end()) {  if ((it->second).size() != 0) { groupFiles = (it->second);  } }
            
                //find group file with the shortest name.  This is because if there is a composite group file it will be the shortest.
                groupFile = groupFiles[0];
                for (int i = 1; i < groupFiles.size(); i++) { if (groupFiles[i].length() < groupFile.length()) { groupFile = groupFiles[i]; } }
            }
            
            inputString = "fasta=" + fastaFile + ", processors=1, name=" + nameFile;
            params->m->mothurOut("\nRunning command: summary.seqs(" + inputString + ")\n");
            
            summarySeqsCommand = new SeqSummaryCommand(inputString);
            summarySeqsCommand->execute();
            
            if (params->m->getControl_pressed()){ break; }
            
            temp = summarySeqsCommand->getOutputFiles();
            mergeOutputFileList(filenames, temp);
            
            delete summarySeqsCommand;
            
            params->m->mothurOut("\n/******************************************/\n");
            
            if (params->append) {
                params->util.appendFiles(fastaFile, params->fasta);
                params->util.appendFiles(nameFile, params->name);
                if (params->makeGroup) { params->util.appendFiles(groupFile, params->group);  }
            }
            
            
            for (it = filenames.begin(); it != filenames.end(); it++) {
                for (int i = 0; i < (it->second).size(); i++) { params->outputNames.push_back((it->second)[i]); params->outputTypes[it->first].push_back((it->second)[i]); }
            }
            
            params->count++;
        }
    }
	catch(exception& e) {
		params->m->errorOut(e, "SffMultipleCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************
long long SffMultipleCommand::createProcesses(vector<string> sffFiles, vector<string> oligosFiles, string fasta, string name, string group){
    try {
#if defined NON_WINDOWS
#else
        //trim.flows, shhh.flows cannot handle multiple processors for windows.
        processors = 1; m->mothurOut("This command can only use 1 processor on Windows platforms, using 1 processors.\n\n");
#endif
        
        current->setMothurCalling(true);

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
		
        //create array of worker threads
        vector<std::thread*> workerThreads;
        vector<sffMultipleData*> data;
        
        //Lauch worker threads
        for (int i = 0; i < processors-1; i++) {
            string extension = toString(i+1);
            
            sffMultipleData* dataBundle = new sffMultipleData(sffFiles, oligosFiles, fasta+extension, name+extension, group+extension, lines[i+1].start, lines[i+1].end);
            
            dataBundle->setVariables(trim, large, allFiles, flip, keepforward, makeGroup, append, largeSize, bdiffs, tdiffs, pdiffs, sdiffs, ldiffs, maxFlows, minFlows, minLength, maxLength, maxHomoP, maxIters, keepFirst, removeLast, signal, noise, cutoff, sigma, flowOrder, lookupFileName, minDelta);
            data.push_back(dataBundle);
            
            workerThreads.push_back(new std::thread(driverSFFMultiple, dataBundle));
        }
        
        sffMultipleData* dataBundle = new sffMultipleData(sffFiles, oligosFiles, fasta, name, group, lines[0].start, lines[0].end);
        dataBundle->setVariables(trim, large, allFiles, flip, keepforward, makeGroup, append, largeSize, bdiffs, tdiffs, pdiffs, sdiffs, ldiffs, maxFlows, minFlows, minLength, maxLength, maxHomoP, maxIters, keepFirst, removeLast, signal, noise, cutoff, sigma, flowOrder, lookupFileName, minDelta);
        
        driverSFFMultiple(dataBundle);
        long long num = dataBundle->count;
        outputNames = dataBundle->outputNames;
        outputTypes = dataBundle->outputTypes;
        
      
        for (int i = 0; i < processors-1; i++) {
            workerThreads[i]->join();
            num += data[i]->count;
            
            outputNames.insert(outputNames.end(), data[i]->outputNames.begin(), data[i]->outputNames.end());
            outputTypes.insert(data[i]->outputTypes.begin(), data[i]->outputTypes.end());
            
            if (append) {
                string extension = toString(i+1);
                util.appendFiles(fasta+extension, fasta);   util.mothurRemove(fasta+extension);
                util.appendFiles(name+extension, name);     util.mothurRemove(name+extension);
                if (makeGroup) { util.appendFiles(group+extension, group);  util.mothurRemove(group+extension); }
            }

            delete data[i];
            delete workerThreads[i];
        }
        delete dataBundle;
        
        current->setMothurCalling(false);
        return num;
    }
	catch(exception& e) {
		m->errorOut(e, "ShhherCommand", "createProcesses");
		exit(1);
	}
}
//**********************************************************************************************************************




