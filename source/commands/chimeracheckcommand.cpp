/*
 *  chimeracheckcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 3/31/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "chimeracheckcommand.h"

//**********************************************************************************************************************
vector<string> ChimeraCheckCommand::setParameters(){	
	try {
		CommandParameter ptemplate("reference", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(ptemplate);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","chimera",false,true,true); parameters.push_back(pfasta);
		CommandParameter pname("name", "InputTypes", "", "", "none", "none", "none","",false,false,true); parameters.push_back(pname);
		CommandParameter psvg("svg", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(psvg);
		CommandParameter pincrement("increment", "Number", "", "10", "", "", "","",false,false); parameters.push_back(pincrement);
		CommandParameter pksize("ksize", "Number", "", "7", "", "", "","",false,false); parameters.push_back(pksize);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);

		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraCheckCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The chimera.check command reads a fastafile and referencefile and outputs potentially chimeric sequences.\n";
		helpString += "This command was created using the algorithms described in CHIMERA_CHECK version 2.7 written by Niels Larsen. \n";
		helpString += "The chimera.check command parameters are fasta, reference, processors, ksize, increment, svg and name.\n";
		helpString += "The fasta parameter allows you to enter the fasta file containing your potentially chimeric sequences, and is required unless you have a valid current fasta file. \n";
		helpString += "The reference parameter allows you to enter a reference file containing known non-chimeric sequences, and is required. \n";
		helpString += "The increment parameter allows you to specify how far you move each window while finding chimeric sequences, default is 10.\n";
		helpString += "The ksize parameter allows you to input kmersize, default is 7. \n";
		helpString += "The svg parameter allows you to specify whether or not you would like a svg file outputted for each query sequence, default is False.\n";
		helpString += "The name parameter allows you to enter a file containing names of sequences you would like .svg files for.\n";
        helpString += "The chimera.check command should be in the following format: \n";
		helpString += "chimera.check(fasta=yourFastaFile, reference=yourTemplateFile, processors=yourProcessors, ksize=yourKmerSize) \n";
		helpString += "Example: chimera.check(fasta=AD.fasta, reference=core_set_aligned,imputed.fasta, processors=4, ksize=8) \n";
			
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraCheckCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "chimera") {  pattern = "[filename],chimeracheck.chimeras"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraCheckCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
ChimeraCheckCommand::ChimeraCheckCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
		vector<string> tempOutNames;
		outputTypes["chimera"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckCommand", "ChimeraCheckCommand");
		exit(1);
	}
}
//***************************************************************************************************************
ChimeraCheckCommand::ChimeraCheckCommand(string option)  {
	try {
		abort = false; calledHelp = false;  
        
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter("chimera.check");
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
			}
			
			vector<string> tempOutNames;
			outputTypes["chimera"] = tempOutNames;
			
            fastafile = validParameter.validFile(parameters, "fasta");
            if (fastafile == "not found") {
                fastafile = current->getFastaFile();
                if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter.\n"); }
                else { 	m->mothurOut("[ERROR]: You have no current fasta file and the fasta parameter is required.\n");  abort = true; }
            }
            else if (fastafile == "not open") { abort = true; }
            else { current->setFastaFile(fastafile); }
            
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";	}
			
            namefile = validParameter.validFile(parameters, "name");
            if (namefile == "not open") { namefile = ""; abort = true; }
            else if (namefile == "not found") {  namefile = "";  }
            else { current->setNameFile(namefile); }
            
			
			//this has to go after save so that if the user sets save=t and provides no reference we abort
			templatefile = validParameter.validFile(parameters, "reference");
			if (templatefile == "not found") {  m->mothurOut("[ERROR]: The reference parameter is a required, aborting.\n"); abort = true;
			}else if (templatefile == "not open") { abort = true; }	
			
			string temp = validParameter.valid(parameters, "ksize");			if (temp == "not found") { temp = "7"; }
			util.mothurConvert(temp, ksize);
			
			temp = validParameter.valid(parameters, "svg");				if (temp == "not found") { temp = "F"; }
			svg = util.isTrue(temp);
			if (namefile != "") { svg = true; }
			
			temp = validParameter.valid(parameters, "increment");		if (temp == "not found") { temp = "10"; }
			util.mothurConvert(temp, increment);			
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckCommand", "ChimeraCheckCommand");
		exit(1);
	}
}
//***************************************************************************************************************

int ChimeraCheckCommand::execute(){
	try{
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
        m->mothurOut("Checking sequences from " + fastafile + " ...\n" );
        
        long start = time(NULL);
        
        numSeqs = checkChimeras();
        
        if (m->getControl_pressed()) { for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	} outputTypes.clear();   return 0; }
       
        m->mothurOut("\nThis method does not determine if a sequence is chimeric, but allows you to make that determination based on the IS values.\n");
        m->mothurOut("\nIt took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences.\n\n");
        
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}	
		m->mothurOutEndLine();
	
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int ChimeraCheckCommand::checkChimeras(){
	try {
        if (outputDir == "") { outputDir = util.hasPath(fastafile);  }//if user entered a file with a path then preserve it
        map<string, string> variables;
        variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(fastafile));
        string outputFileName = getOutputFileName("chimera", variables);
        outputNames.push_back(outputFileName); outputTypes["chimera"].push_back(outputFileName);
        
        MothurChimera* chimera = new ChimeraCheckRDP(fastafile, templatefile, namefile, svg, increment, ksize, outputDir);
        
		ofstream out;
		util.openOutputFile(outputFileName, out);
		
		ofstream out2;
		
		ifstream inFASTA;
		util.openInputFile(fastafile, inFASTA);

		int count = 0;
	
		while (!inFASTA.eof()) {

            if (m->getControl_pressed()) {	break;	}
		
			Sequence* candidateSeq = new Sequence(inFASTA);  util.gobble(inFASTA);
				
			if (candidateSeq->getName() != "") { //incase there is a commented sequence at the end of a file
				//find chimeras
				chimera->getChimeras(candidateSeq);
				
				if (m->getControl_pressed()) {	delete candidateSeq; return 1;	}
	
				//print results
				chimera->print(out, out2);
                count++;
			}
			delete candidateSeq;
			
			//report progress
			if((count) % 100 == 0){	m->mothurOutJustToScreen("Processing sequence: " + toString(count) + "\n");		}
		}
		//report progress
		if((count) % 100 != 0){	m->mothurOutJustToScreen("Processing sequence: " + toString(count) + "\n");	}
		
		out.close();
		inFASTA.close();
        
        delete chimera;
				
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCheckCommand", "checkChimeras");
		exit(1);
	}
}
/**************************************************************************************************/


