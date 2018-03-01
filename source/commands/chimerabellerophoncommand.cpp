/*
 *  chimerabellerophoncommand.cpp
 *  Mothur
 *
 *  Created by westcott on 4/1/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "chimerabellerophoncommand.h"
#include "bellerophon.h"

//**********************************************************************************************************************
vector<string> ChimeraBellerophonCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none","none","none","chimera-accnos",false,true,true); parameters.push_back(pfasta);
		CommandParameter pfilter("filter", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pfilter);
		CommandParameter pcorrection("correction", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pcorrection);
		CommandParameter pwindow("window", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pwindow);
		CommandParameter pincrement("increment", "Number", "", "25", "", "", "","",false,false); parameters.push_back(pincrement);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
		
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraBellerophonCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraBellerophonCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The chimera.bellerophon command reads a fastafile and creates list of potentially chimeric sequences.\n";
		helpString += "The chimera.bellerophon command parameters are fasta, filter, correction, processors, window, increment. The fasta parameter is required, unless you have a valid current file.\n";
		helpString += "The fasta parameter is required.  You may enter multiple fasta files by separating their names with dashes. ie. fasta=abrecovery.fasta-amzon.fasta \n";
		helpString += "The filter parameter allows you to specify if you would like to apply a vertical and 50% soft filter, default=false. \n";
		helpString += "The correction parameter allows you to put more emphasis on the distance between highly similar sequences and less emphasis on the differences between remote homologs, default=true.\n";
		helpString += "The window parameter allows you to specify the window size for searching for chimeras, default is 1/4 sequence length. \n";
		helpString += "The increment parameter allows you to specify how far you move each window while finding chimeric sequences, default is 25.\n";
		helpString += "chimera.bellerophon(fasta=yourFastaFile, filter=yourFilter, correction=yourCorrection, processors=yourProcessors) \n";
		helpString += "Example: chimera.bellerophon(fasta=AD.align, filter=True, correction=true, window=200) \n";
		helpString += "Note: No spaces between parameter labels (i.e. fasta), '=' and parameters (i.e.yourFastaFile).\n";	
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraBellerophonCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraBellerophonCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "chimera") {  pattern = "[filename],bellerophon.chimeras"; } 
        else if (type == "accnos") {  pattern = "[filename],bellerophon.accnos"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraBellerophonCommand", "getOutputPattern");
        exit(1);
    }
}
//**********************************************************************************************************************
ChimeraBellerophonCommand::ChimeraBellerophonCommand(){	
	try {
		abort = true; calledHelp = true; 
		setParameters();
		vector<string> tempOutNames;
		outputTypes["chimera"] = tempOutNames;
		outputTypes["accnos"] = tempOutNames;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraBellerophonCommand", "ChimeraBellerophonCommand");
		exit(1);
	}
}
//***************************************************************************************************************
ChimeraBellerophonCommand::ChimeraBellerophonCommand(string option)  {
	try {
		abort = false; calledHelp = false;   
		
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
		
		else {
			
			vector<string> myArray = setParameters();
			
			OptionParser parser(option);
			map<string,string> parameters = parser.getParameters();
			
			ValidParameters validParameter("chimera.bellerophon");
			map<string,string>::iterator it;
			
			//check to make sure all parameters are valid for command
			for (it = parameters.begin(); it != parameters.end(); it++) { 
				if (validParameter.isValidParameter(it->first, myArray, it->second) != true) {  abort = true;  }
			}
			
			//initialize outputTypes
			vector<string> tempOutNames;
			outputTypes["chimera"] = tempOutNames;
			outputTypes["accnos"] = tempOutNames;
		
			//if the user changes the input directory command factory will send this info to us in the output parameter 
			string inputDir = validParameter.valid(parameters, "inputdir");
			if (inputDir == "not found"){	inputDir = "";		}
			
			fastafile = validParameter.valid(parameters, "fasta");
			if (fastafile == "not found") { 				
				//if there is a current fasta file, use it
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

			string temp;
			temp = validParameter.valid(parameters, "filter");			if (temp == "not found") { temp = "F"; }
			filter = util.isTrue(temp);
			
			temp = validParameter.valid(parameters, "correction");		if (temp == "not found") { temp = "T"; }
			correction = util.isTrue(temp);
			
			temp = validParameter.valid(parameters, "window");			if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, window);
			
			temp = validParameter.valid(parameters, "increment");		if (temp == "not found") { temp = "25"; }
			util.mothurConvert(temp, increment);
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraBellerophonCommand", "ChimeraBellerophonCommand");
		exit(1);
	}
}
//***************************************************************************************************************
int ChimeraBellerophonCommand::execute(){
	try{
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
		for (int i = 0; i < fastaFileNames.size(); i++) {
			
			m->mothurOut("Checking sequences from " + fastaFileNames[i] + " ..." ); m->mothurOutEndLine();
			
			long start = time(NULL);	
			
			chimera = new Bellerophon(fastaFileNames[i], filter, correction, window, increment, outputDir);	
			
			if (outputDir == "") { outputDir = util.hasPath(fastaFileNames[i]);  }//if user entered a file with a path then preserve it	
            map<string, string> variables; 
            variables["[filename]"] = outputDir + util.getRootName(util.getSimpleName(fastaFileNames[i]));
			string outputFileName = getOutputFileName("chimera", variables);
			string accnosFileName = getOutputFileName("accnos", variables);
			
			chimera->getChimeras();
			
			if (m->getControl_pressed()) { delete chimera; for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);	} outputTypes.clear(); return 0;	}
					
			ofstream out;
			util.openOutputFile(outputFileName, out);
			
			ofstream out2;
			util.openOutputFile(accnosFileName, out2);
			
			numSeqs = chimera->print(out, out2, "");
			out.close();
			out2.close(); 
						
			if (m->getControl_pressed()) { util.mothurRemove(accnosFileName); util.mothurRemove(outputFileName); for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);	} outputTypes.clear(); delete chimera;	return 0;	}
			
			m->mothurOutEndLine(); m->mothurOut("It took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences.");	m->mothurOutEndLine(); m->mothurOutEndLine();
			
			outputNames.push_back(outputFileName);  outputTypes["chimera"].push_back(outputFileName);
			outputNames.push_back(accnosFileName);  outputTypes["accnos"].push_back(accnosFileName);
			
			delete chimera;
		}
		
		//set accnos file as new current accnosfile
		string currentName = "";
		itTypes = outputTypes.find("accnos");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setAccnosFile(currentName); }
		}
		
		m->mothurOut("\nOutput File Names: \n"); 
		for (int i = 0; i < outputNames.size(); i++) {	m->mothurOut(outputNames[i] +"\n"); 	} m->mothurOutEndLine();
		
		return 0;
				
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraBellerophonCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

