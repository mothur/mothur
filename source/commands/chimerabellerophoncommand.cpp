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
#include "removeseqscommand.h"


//**********************************************************************************************************************
vector<string> ChimeraBellerophonCommand::setParameters(){	
	try {
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none","none","none","chimera-accnos",false,true,true); parameters.push_back(pfasta);
		CommandParameter pfilter("filter", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pfilter);
		CommandParameter pcorrection("correction", "Boolean", "", "T", "", "", "","",false,false); parameters.push_back(pcorrection);
		CommandParameter pwindow("window", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pwindow);
		CommandParameter pincrement("increment", "Number", "", "25", "", "", "","",false,false); parameters.push_back(pincrement);
        CommandParameter premovechimeras("removechimeras", "Boolean", "", "t", "", "", "","fasta",false,false); parameters.push_back(premovechimeras);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);
        
        abort = false; calledHelp = false;
       
        vector<string> tempOutNames;
        outputTypes["chimera"] = tempOutNames;
        outputTypes["accnos"] = tempOutNames;
        outputTypes["fasta"] = tempOutNames;

		
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
		helpString += "The filter parameter allows you to specify if you would like to apply a vertical and 50% soft filter, default=false. \n";
		helpString += "The correction parameter allows you to put more emphasis on the distance between highly similar sequences and less emphasis on the differences between remote homologs, default=true.\n";
		helpString += "The window parameter allows you to specify the window size for searching for chimeras, default is 1/4 sequence length. \n";
		helpString += "The increment parameter allows you to specify how far you move each window while finding chimeric sequences, default is 25.\n";
        helpString += "The removechimeras parameter allows you to indicate you would like to automatically remove the sequences that are flagged as chimeric. Default=t.\n";

		helpString += "chimera.bellerophon(fasta=yourFastaFile, filter=yourFilter, correction=yourCorrection, processors=yourProcessors) \n";
		helpString += "Example: chimera.bellerophon(fasta=AD.align, filter=True, correction=true, window=200) \n";
			
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
        else if (type == "fasta")   {  pattern = "[filename],bellerophon.fasta";     }
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraBellerophonCommand", "getOutputPattern");
        exit(1);
    }
}
//***************************************************************************************************************
ChimeraBellerophonCommand::ChimeraBellerophonCommand(string option) : Command()  {
	try {
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();
            
            ValidParameters validParameter;
            fastafile = validParameter.validFile(parameters, "fasta");
            if (fastafile == "not found") {
                fastafile = current->getFastaFile();
                if (fastafile != "") { m->mothurOut("Using " + fastafile + " as input file for the fasta parameter.\n"); }
                else { 	m->mothurOut("[ERROR]: You have no current fasta file and the fasta parameter is required.\n");  abort = true; }
            }
            else if (fastafile == "not open") { abort = true; }
            else { current->setFastaFile(fastafile); }
			
			string temp;
			temp = validParameter.valid(parameters, "filter");			if (temp == "not found") { temp = "F"; }
			filter = util.isTrue(temp);
			
			temp = validParameter.valid(parameters, "correction");		if (temp == "not found") { temp = "T"; }
			correction = util.isTrue(temp);
			
			temp = validParameter.valid(parameters, "window");			if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, window);
			
			temp = validParameter.valid(parameters, "increment");		if (temp == "not found") { temp = "25"; }
			util.mothurConvert(temp, increment);
            
            temp = validParameter.valid(parameters, "removechimeras");            if (temp == "not found") { temp = "t"; }
            removeChimeras = util.isTrue(temp);
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
		
        m->mothurOut("Checking sequences from " + fastafile + " ...\n" );
        
        long start = time(nullptr);

        MothurChimera* chimera = new Bellerophon(fastafile, filter, correction, window, increment, outputdir);
        
        chimera->getChimeras();
        
        if (m->getControl_pressed()) { delete chimera;  return 0;	}
        
        if (outputdir == "") { outputdir = util.hasPath(fastafile);  }
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(fastafile));
        string outputFileName = getOutputFileName("chimera", variables);
        string accnosFileName = getOutputFileName("accnos", variables);
        
        ofstream out; util.openOutputFile(outputFileName, out); outputNames.push_back(outputFileName);  outputTypes["chimera"].push_back(outputFileName);
        ofstream out2; util.openOutputFile(accnosFileName, out2); outputNames.push_back(accnosFileName);  outputTypes["accnos"].push_back(accnosFileName);
        
        //print results
        numSeqs = chimera->print(out, out2, "");
        
        out.close(); out2.close();
        
        delete chimera;
        
        if (m->getControl_pressed()) {  for (int i = 0; i < outputNames.size(); i++) {	util.mothurRemove(outputNames[i]);	} outputTypes.clear(); 	return 0;	}
        
        m->mothurOut("\nIt took " + toString(time(nullptr) - start) + " secs to check " + toString(numSeqs) + " sequences.\n\n");
        
        if (removeChimeras) {
            if (!util.isBlank(accnosFileName)) {
                m->mothurOut("\nRemoving chimeras from your input files:\n");
                
                string inputString = "fasta=" + fastafile + ", accnos=" + accnosFileName;
                
                m->mothurOut("/******************************************/\n");
                m->mothurOut("Running command: remove.seqs(" + inputString + ")\n");
                current->setMothurCalling(true);
                
                Command* removeCommand = new RemoveSeqsCommand(inputString);
                removeCommand->execute();
                
                map<string, vector<string> > filenames = removeCommand->getOutputFiles();
                
                delete removeCommand;
                current->setMothurCalling(false);
                m->mothurOut("/******************************************/\n");
                
                map<string, string> variables;
                variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(fastafile));
                string currentName = getOutputFileName("fasta", variables);

                util.renameFile(filenames["fasta"][0], currentName);
                util.mothurRemove(filenames["fasta"][0]);
                
                outputNames.push_back(currentName); outputTypes["fasta"].push_back(currentName);
            }else { m->mothurOut("\nNo chimeras found, skipping remove.seqs.\n"); }
        }

		//set accnos file as new current accnosfile
		string currentName = "";
		itTypes = outputTypes.find("accnos");
		if (itTypes != outputTypes.end()) {
			if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setAccnosFile(currentName); }
		}
		
        itTypes = outputTypes.find("fasta");
        if (itTypes != outputTypes.end()) {
            if ((itTypes->second).size() != 0) { currentName = (itTypes->second)[0]; current->setFastaFile(currentName); }
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

