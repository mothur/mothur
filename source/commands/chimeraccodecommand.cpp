/*
 *  chimeraccodecommand.cpp
 *  Mothur
 *
 *  Created by westcott on 3/30/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "chimeraccodecommand.h"
#include "ccode.h"

//**********************************************************************************************************************
vector<string> ChimeraCcodeCommand::setParameters(){	
	try {
		CommandParameter ptemplate("reference", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(ptemplate);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","chimera-mapinfo-accnos",false,true,true); parameters.push_back(pfasta);
		CommandParameter pfilter("filter", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pfilter);
		CommandParameter pwindow("window", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pwindow);
		CommandParameter pnumwanted("numwanted", "Number", "", "20", "", "", "","",false,false); parameters.push_back(pnumwanted);
		CommandParameter pmask("mask", "String", "", "", "", "", "","",false,false); parameters.push_back(pmask);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCcodeCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraCcodeCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The chimera.ccode command reads a fastafile and referencefile and outputs potentially chimeric sequences.\n";
		helpString += "This command was created using the algorithms described in the 'Evaluating putative chimeric sequences from PCR-amplified products' paper by Juan M. Gonzalez, Johannes Zimmerman and Cesareo Saiz-Jimenez.\n";
		helpString += "The chimera.ccode command parameters are fasta, reference, filter, mask, processors, window and numwanted.\n";
		helpString += "The fasta parameter allows you to enter the fasta file containing your potentially chimeric sequences, and is required unless you have a valid current fasta file. \n";
		helpString += "You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amzon.fasta \n";
		helpString += "The reference parameter allows you to enter a reference file containing known non-chimeric sequences, and is required. \n";
		helpString += "The filter parameter allows you to specify if you would like to apply a vertical and 50% soft filter. \n";
		helpString += "The mask parameter allows you to specify a file containing one sequence you wish to use as a mask for the your sequences. \n";
		helpString += "The window parameter allows you to specify the window size for searching for chimeras. \n";
		helpString += "The numwanted parameter allows you to specify how many sequences you would each query sequence compared with.\n";
		helpString += "The chimera.ccode command should be in the following format: \n";
		helpString += "chimera.ccode(fasta=yourFastaFile, reference=yourTemplate) \n";
		helpString += "Example: chimera.ccode(fasta=AD.align, reference=core_set_aligned.imputed.fasta) \n";
			
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCcodeCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraCcodeCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "chimera") {  pattern = "[filename],[tag],ccode.chimeras-[filename],ccode.chimeras"; } 
        else if (type == "accnos") {  pattern = "[filename],[tag],ccode.accnos-[filename],ccode.accnos"; } 
        else if (type == "mapinfo") {  pattern =  "[filename],mapinfo"; }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraCcodeCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
ChimeraCcodeCommand::ChimeraCcodeCommand(){	
	try {
		abort = true; calledHelp = true;
		setParameters();
		vector<string> tempOutNames;
		outputTypes["chimera"] = tempOutNames;
		outputTypes["mapinfo"] = tempOutNames;
		outputTypes["accnos"] = tempOutNames;
        
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCcodeCommand", "ChimeraCcodeCommand");
		exit(1);
	}
}
//***************************************************************************************************************
ChimeraCcodeCommand::ChimeraCcodeCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter("chimera.ccode");
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (!validParameter.isValidParameter(it->first, myArray, it->second)) {  abort = true;  }
			}
			
			vector<string> tempOutNames;
			outputTypes["chimera"] = tempOutNames;
			outputTypes["mapinfo"] = tempOutNames;
			outputTypes["accnos"] = tempOutNames;
            
			
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");
			if (inputDir == "not found"){	inputDir = "";		}
			else {
				string path;
				it = parameters.find("reference");
				//user has given a template file
				if(it != parameters.end()){ 
					path = util.hasPath(it->second);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	parameters["reference"] = inputDir + it->second;		}
				}
			}

			//check for required parameters
			fastafile = validParameter.valid(parameters, "fasta");
			if (fastafile == "not found") { 				//if there is a current fasta file, use it
				string filename = current->getFastaFile(); 
				if (filename != "") { fastaFileNames.push_back(filename); m->mothurOut("Using " + filename + " as input file for the fasta parameter."); m->mothurOutEndLine(); }
				else { 	m->mothurOut("You have no current fastafile and the fasta parameter is required."); m->mothurOutEndLine(); abort = true; }
			}else { 
				util.splitAtDash(fastafile, fastaFileNames);
				
				//go through files and make sure they are good, if not, then disregard them
				for (int i = 0; i < fastaFileNames.size(); i++) {
					
					bool ignore = false;
					if (fastaFileNames[i] == "current") { 
						fastaFileNames[i] = current->getFastaFile(); 
						if (fastaFileNames[i] != "") {  m->mothurOut("Using " + fastaFileNames[i] + " as input file for the fasta parameter where you had given current."); m->mothurOutEndLine(); }
						else { 	
							m->mothurOut("You have no current fastafile, ignoring current."); m->mothurOutEndLine(); ignore=true; 
							//erase from file list
							fastaFileNames.erase(fastaFileNames.begin()+i);
							i--;
						}
					}
					
                    if (!ignore) {
                        if (util.checkLocations(fastaFileNames[i], current->getLocations())) { current->setFastaFile(fastaFileNames[i]); }
                        else { fastaFileNames.erase(fastaFileNames.begin()+i); i--; } //erase from file list
                    }
                }
				
				//make sure there is at least one valid file left
				if (fastaFileNames.size() == 0) { m->mothurOut("no valid files."); m->mothurOutEndLine(); abort = true; }
			}
			
			//if the user changes the output directory command factory will send this info to us in the output parameter 
			outputDir = validParameter.valid(parameters, "outputdir");		if (outputDir == "not found"){	outputDir = "";	}
			
			maskfile = validParameter.valid(parameters, "mask");
			if (maskfile == "not found") { maskfile = "";  }	
			else if (maskfile != "default")  { 
				if (inputDir != "") {
					string path = util.hasPath(maskfile);
					//if the user has not given a path then, add inputdir. else leave path alone.
					if (path == "") {	maskfile = inputDir + maskfile;		}
				}

				ifstream in;
				bool ableToOpen = util.openInputFile(maskfile, in);
				if (!ableToOpen) { abort = true; }
				in.close();
			}
			
			string temp;
			temp = validParameter.valid(parameters, "filter");			if (temp == "not found") { temp = "F"; }
			filter = util.isTrue(temp);
			
			temp = validParameter.valid(parameters, "window");			if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, window);
			
			temp = validParameter.valid(parameters, "numwanted");		if (temp == "not found") { temp = "20"; }
			util.mothurConvert(temp, numwanted);
			
			//this has to go after save so that if the user sets save=t and provides no reference we abort
			templatefile = validParameter.validFile(parameters, "reference");
			if (templatefile == "not found") { m->mothurOut("[ERROR]: The reference parameter is a required, aborting.\n"); abort = true;
			}else if (templatefile == "not open") { abort = true; }

		}
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCcodeCommand", "ChimeraCcodeCommand");
		exit(1);
	}
}
//***************************************************************************************************************
int ChimeraCcodeCommand::execute(){
	try{
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		for (int s = 0; s < fastaFileNames.size(); s++) {
				
			m->mothurOut("Checking sequences from " + fastaFileNames[s] + " ..." ); m->mothurOutEndLine();
		
			long start = time(NULL);	
			
			//set user options
			if (maskfile == "default") { m->mothurOut("I am using the default 236627 EU009184.1 Shigella dysenteriae str. FBD013."); m->mothurOutEndLine();  }

			chimera = new Ccode(fastaFileNames[s], templatefile, filter, maskfile, window, numwanted, outputDir);	
			
			//is your template aligned?
			if (chimera->getUnaligned()) { m->mothurOut("Your template sequences are different lengths, please correct."); m->mothurOutEndLine(); delete chimera; return 0; }
			templateSeqsLength = chimera->getLength();
			
			if (outputDir == "") { outputDir = util.hasPath(fastaFileNames[s]);  }//if user entered a file with a path then preserve it
			string outputFileName, accnosFileName;
            map<string, string> variables; 
            variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(fastaFileNames[s]));
            string mapInfo = getOutputFileName("mapinfo", variables);
			if (maskfile != "") { variables["[tag]"] = maskfile; }
            outputFileName = getOutputFileName("chimera", variables);
            accnosFileName = getOutputFileName("accnos", variables);
						
			if (m->getControl_pressed()) { delete chimera;  for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	} outputTypes.clear(); return 0;	}
			
			ofstream outHeader;
			string tempHeader = outputDir + util.getRootName(util.getSimpleName(fastaFileNames[s])) + maskfile + "ccode.chimeras.tempHeader";
			util.openOutputFile(tempHeader, outHeader);
			
			outHeader << "For full window mapping info refer to " << mapInfo << endl << endl; outHeader.close();
							
            numSeqs = driver(outputFileName, fastaFileNames[s], accnosFileName);
					
            if (m->getControl_pressed()) { util.mothurRemove(outputFileName); util.mothurRemove(tempHeader); util.mothurRemove(accnosFileName); for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  outputTypes.clear();  delete chimera; return 0; }

            util.appendFiles(outputFileName, tempHeader);
		
			util.mothurRemove(outputFileName);
			rename(tempHeader.c_str(), outputFileName.c_str());
		
			delete chimera;
			
			outputNames.push_back(outputFileName); outputTypes["chimera"].push_back(outputFileName);
			outputNames.push_back(mapInfo);	outputTypes["mapinfo"].push_back(mapInfo);
			outputNames.push_back(accnosFileName); outputTypes["accnos"].push_back(accnosFileName);
			 
						
			m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences.");	m->mothurOutEndLine();
		}
		
		
		//set accnos file as new current accnosfile
		string currentName = "";
		itTypes = outputTypes.find("accnos");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setAccnosFile(currentName); }
		}
		
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i]); m->mothurOutEndLine();	}	
		m->mothurOutEndLine();
		
		return 0;
		
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCcodeCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int ChimeraCcodeCommand::driver(string outputFName, string filename, string accnos){
	try {
		ofstream out;
		util.openOutputFile(outputFName, out);
		
		ofstream out2;
		util.openOutputFile(accnos, out2);
		
		ifstream inFASTA;
		util.openInputFile(filename, inFASTA);

		int count = 0;
	
		while (!inFASTA.eof()) {
		
			if (m->getControl_pressed()) {	return 1;	}
		
			Sequence* candidateSeq = new Sequence(inFASTA);  util.gobble(inFASTA);
				
			if (candidateSeq->getName() != "") { //incase there is a commented sequence at the end of a file
				
				if (candidateSeq->getAligned().length() != templateSeqsLength) {  
					m->mothurOut(candidateSeq->getName() + " is not the same length as the template sequences. Skipping."); m->mothurOutEndLine();
				}else{
					//find chimeras
					chimera->getChimeras(candidateSeq);
					
					if (m->getControl_pressed()) {	delete candidateSeq; return 1;	}
		
					//print results
					chimera->print(out, out2);
				}
				count++;
			}
			delete candidateSeq;
			
						
			//report progress
			if((count) % 100 == 0){	m->mothurOutJustToScreen("Processing sequence: " + toString(count) + "\n");		}
		}
		//report progress
		if((count) % 100 != 0){	m->mothurOutJustToScreen("Processing sequence: " + toString(count) + "\n");	}
		
		out.close();
		out2.close();
		inFASTA.close();
				
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraCcodeCommand", "driver");
		exit(1);
	}
}
//**********************************************************************************************************************

