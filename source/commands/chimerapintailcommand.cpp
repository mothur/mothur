/*
 *  chimerapintailcommand.cpp
 *  Mothur
 *
 *  Created by westcott on 4/1/10.
 *  Copyright 2010 Schloss Lab. All rights reserved.
 *
 */

#include "chimerapintailcommand.h"
#include "pintail.h"


//**********************************************************************************************************************
vector<string> ChimeraPintailCommand::setParameters(){	
	try {
		CommandParameter ptemplate("reference", "InputTypes", "", "", "none", "none", "none","",false,true,true); parameters.push_back(ptemplate);
		CommandParameter pfasta("fasta", "InputTypes", "", "", "none", "none", "none","chimera-accnos",false,true,true); parameters.push_back(pfasta);
		CommandParameter pconservation("conservation", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pconservation);
		CommandParameter pquantile("quantile", "InputTypes", "", "", "none", "none", "none","",false,false); parameters.push_back(pquantile);
		CommandParameter pfilter("filter", "Boolean", "", "F", "", "", "","",false,false); parameters.push_back(pfilter);
		CommandParameter pwindow("window", "Number", "", "0", "", "", "","","",false,false); parameters.push_back(pwindow);
		CommandParameter pincrement("increment", "Number", "", "25", "", "", "","",false,false); parameters.push_back(pincrement);
		CommandParameter pmask("mask", "String", "", "", "", "", "","",false,false); parameters.push_back(pmask);
		CommandParameter pseed("seed", "Number", "", "0", "", "", "","",false,false); parameters.push_back(pseed);
        CommandParameter pinputdir("inputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(pinputdir);
		CommandParameter poutputdir("outputdir", "String", "", "", "", "", "","",false,false); parameters.push_back(poutputdir);

        abort = false; calledHelp = false;
        
        vector<string> tempOutNames;
        outputTypes["chimera"] = tempOutNames;
        outputTypes["accnos"] = tempOutNames;
        
		vector<string> myArray;
		for (int i = 0; i < parameters.size(); i++) {	myArray.push_back(parameters[i].name);		}
		return myArray;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPintailCommand", "setParameters");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraPintailCommand::getHelpString(){	
	try {
		string helpString = "";
		helpString += "The chimera.pintail command reads a fastafile and referencefile and outputs potentially chimeric sequences.\n";
		helpString += "This command was created using the algorithms described in the 'At Least 1 in 20 16S rRNA Sequence Records Currently Held in the Public Repositories is Estimated To Contain Substantial Anomalies' paper by Kevin E. Ashelford 1, Nadia A. Chuzhanova 3, John C. Fry 1, Antonia J. Jones 2 and Andrew J. Weightman 1.\n";
		helpString += "The chimera.pintail command parameters are fasta, reference, filter, mask, processors, window, increment, conservation and quantile.\n";
		helpString += "The fasta parameter allows you to enter the fasta file containing your potentially chimeric sequences, and is required unless you have a valid current fasta file. \n";
		helpString += "The reference parameter allows you to enter a reference file containing known non-chimeric sequences, and is required. \n";
		helpString += "The filter parameter allows you to specify if you would like to apply a vertical and 50% soft filter. \n";
		helpString += "The mask parameter allows you to specify a file containing one sequence you wish to use as a mask for the your sequences, by default no mask is applied.  You can apply an ecoli mask by typing, mask=default. \n";
		helpString += "The window parameter allows you to specify the window size for searching for chimeras, default=300. \n";
		helpString += "The increment parameter allows you to specify how far you move each window while finding chimeric sequences, default=25.\n";
		helpString += "The conservation parameter allows you to enter a frequency file containing the highest bases frequency at each place in the alignment.\n";
		helpString += "The quantile parameter allows you to enter a file containing quantiles for a template files sequences, if you use the filter the quantile file generated becomes unique to the fasta file you used.\n";
		helpString += "The chimera.pintail command should be in the following format: \n";
		helpString += "chimera.pintail(fasta=yourFastaFile, reference=yourTemplate) \n";
		helpString += "Example: chimera.pintail(fasta=AD.align, reference=silva.bacteria.fasta) \n";
			
		return helpString;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPintailCommand", "getHelpString");
		exit(1);
	}
}
//**********************************************************************************************************************
string ChimeraPintailCommand::getOutputPattern(string type) {
    try {
        string pattern = "";
        
        if (type == "chimera") {  pattern = "[filename],[tag],pintail.chimeras-[filename],pintail.chimeras"; } 
        else if (type == "accnos") {  pattern = "[filename],[tag],pintail.accnos-[filename],pintail.accnos"; } 
        else { m->mothurOut("[ERROR]: No definition for type " + type + " output pattern.\n"); m->setControl_pressed(true);  }
        
        return pattern;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraPintailCommand", "getOutputPattern");
        exit(1);
    }
}
//***************************************************************************************************************
ChimeraPintailCommand::ChimeraPintailCommand(string option)  {
	try {
		//allow user to run help
		if(option == "help") { help(); abort = true; calledHelp = true; }
		else if(option == "citation") { citation(); abort = true; calledHelp = true;}
        else if(option == "category") {  abort = true; calledHelp = true;  }
		
		else {
			OptionParser parser(option, setParameters());
			map<string,string> parameters = parser.getParameters();

			//check for required parameters
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
						
			temp = validParameter.valid(parameters, "window");			if (temp == "not found") { temp = "0"; }
			util.mothurConvert(temp, window);
			
			temp = validParameter.valid(parameters, "increment");		if (temp == "not found") { temp = "25"; }
			util.mothurConvert(temp, increment);
			
			
			//this has to go after save so that if the user sets save=t and provides no reference we abort
			templatefile = validParameter.validFile(parameters, "reference");
			if (templatefile == "not found") { m->mothurOut("[ERROR]: The reference parameter is a required, aborting.\n"); abort = true;
			}else if (templatefile == "not open") { abort = true; }
			
			
			maskfile = validParameter.validPath(parameters, "mask");
			if (maskfile == "not found") { maskfile = "";  }	
			else if (maskfile != "default")  {
                if (util.checkLocations(maskfile, current->getLocations())) {  }
                else { m->mothurOut("Unable to open " + maskfile + ".\n"); abort = true; } //erase from file list
            }else if (maskfile == "default") { m->mothurOut("Using the default 236627 EU009184.1 Shigella dysenteriae str. FBD013.\n");  }

			 
			
            if (outputdir == "") { outputdir = util.hasPath(fastafile);  }
		
			consfile = validParameter.validFile(parameters, "conservation");
			if (consfile == "not open") { abort = true; }
			else if (consfile == "not found") { 
				consfile = "";  
				//check for consfile
				string tempConsFile = util.getRootName(inputDir + util.getSimpleName(templatefile)) + "freq";
				ifstream FileTest(tempConsFile.c_str());
				if(FileTest){	
                    string line = util.getline(FileTest);
                    bool GoodFile = util.checkReleaseVersion(line, current->getVersion());
					if (GoodFile) {  
						m->mothurOut("I found " + tempConsFile + " in your input file directory. I will use it to save time.\n"); consfile = tempConsFile;  FileTest.close();
					}
				}else {
					string tempConsFile = current->getDefaultPath() + util.getRootName(util.getSimpleName(templatefile)) + "freq";
					ifstream FileTest2(tempConsFile.c_str());
					if(FileTest2){
                        string line = util.getline(FileTest2);
						bool GoodFile = util.checkReleaseVersion(line, current->getVersion());
						if (GoodFile) {  
							m->mothurOut("I found " + tempConsFile + " in your input file directory. I will use it to save time.\n"); consfile = tempConsFile;  FileTest2.close();
						}
					}
				}
			}	
			
			quanfile = validParameter.validFile(parameters, "quantile");
			if (quanfile == "not open") { abort = true; }
			else if (quanfile == "not found") { quanfile = ""; }
		}
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPintailCommand", "ChimeraPintailCommand");
		exit(1);
	}
}
//***************************************************************************************************************

int ChimeraPintailCommand::execute(){
	try{
		
		if (abort) { if (calledHelp) { return 0; }  return 2;	}
		
        m->mothurOut("Checking sequences from " + fastafile + " ...\n" );
        
        long start = time(NULL);
        
        //check for quantile to save the time
        lookForShortcutFiles(templatefile);
        
        numSeqs = checkChiemras();
        
        if (m->getControl_pressed()) { outputTypes.clear();  for (int j = 0; j < outputNames.size(); j++) {	util.mothurRemove(outputNames[j]);	}  return 0; }
        
        m->mothurOut("\n\nIt took " + toString(time(NULL) - start) + " secs to check " + toString(numSeqs) + " sequences.\n");
        
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
		m->errorOut(e, "ChimeraPintailCommand", "execute");
		exit(1);
	}
}
//**********************************************************************************************************************

int ChimeraPintailCommand::checkChiemras(){
	try {
        MothurChimera* chimera = new Pintail(fastafile, templatefile, filter, maskfile, consfile, quanfile, window, increment, outputdir, current->getVersion());
        
        if (m->getControl_pressed()) { delete chimera;  return 0;	}
        
        if (chimera->getUnaligned()) { m->mothurOut("[ERROR]: Your reference sequences are unaligned, please correct.\n"); delete chimera; return 0; }
        
        templateSeqsLength = chimera->getLength();
        
        string outputFileName, accnosFileName;
        map<string, string> variables;
        variables["[filename]"] = outputdir + util.getRootName(util.getSimpleName(fastafile));
        if (maskfile != "") { variables["[tag]"] = util.getSimpleName(util.getRootName(maskfile)); }
        outputFileName = getOutputFileName("chimera", variables);
        accnosFileName = getOutputFileName("accnos", variables);

        outputNames.push_back(outputFileName); outputTypes["chimera"].push_back(outputFileName);
        outputNames.push_back(accnosFileName); outputTypes["accnos"].push_back(accnosFileName);

		ofstream out;
		util.openOutputFile(outputFileName, out);
		
		ofstream out2;
		util.openOutputFile(accnosFileName, out2);
		
		ifstream inFASTA;
		util.openInputFile(fastafile, inFASTA);

		int count = 0;
	
		while (!inFASTA.eof()) {
				
			if (m->getControl_pressed()) {	break;	}
		
			Sequence* candidateSeq = new Sequence(inFASTA);  util.gobble(inFASTA);
				
			if (candidateSeq->getName() != "") { //incase there is a commented sequence at the end of a file
				
				if (candidateSeq->getAligned().length() != templateSeqsLength)  {  //chimeracheck does not require seqs to be aligned
					m->mothurOut("[WARNING]: " + candidateSeq->getName() + " is not the same length as the template sequences. Skipping.\n");
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
			if((count) % 100 == 0){	m->mothurOutJustToScreen("Processing sequence: " + toString(count) + "\n"); 		}
		}
		//report progress
		if((count) % 100 != 0){	m->mothurOutJustToScreen("Processing sequence: " + toString(count) + "\n"); 		}
		
		out.close();
		out2.close();
		inFASTA.close();
        delete chimera;
				
		return count;
	}
	catch(exception& e) {
		m->errorOut(e, "ChimeraPintailCommand", "checkChiemras");
		exit(1);
	}
}
//**********************************************************************************************************************

int ChimeraPintailCommand::lookForShortcutFiles(string baseName){
    try {
        string tempQuan = "";
        if ((!filter) && (maskfile == "")) {
            tempQuan = inputDir + util.getRootName(util.getSimpleName(baseName)) + "pintail.quan";
        }else if ((!filter) && (maskfile != "")) {
            tempQuan = inputDir + util.getRootName(util.getSimpleName(baseName)) + "pintail.masked.quan";
        }else if ((filter) && (maskfile != "")) {
            tempQuan = inputDir + util.getRootName(util.getSimpleName(baseName)) + "pintail.filtered." + util.getSimpleName(util.getRootName(fastafile)) + "masked.quan";
        }else if ((filter) && (maskfile == "")) {
            tempQuan = inputDir + util.getRootName(util.getSimpleName(baseName)) + "pintail.filtered." + util.getSimpleName(util.getRootName(fastafile)) + "quan";
        }
        
        ifstream FileTest(tempQuan.c_str());
        if(FileTest){
            string line = util.getline(FileTest);
            bool GoodFile = util.checkReleaseVersion(line, current->getVersion());
            if (GoodFile) {
                m->mothurOut("I found " + tempQuan + " in your input file directory. I will use it to save time.\n"); quanfile = tempQuan;  FileTest.close();
            }
        }else {
            string tryPath = current->getDefaultPath();
            string tempQuan = "";
            if ((!filter) && (maskfile == "")) {
                tempQuan = tryPath + util.getRootName(util.getSimpleName(baseName)) + "pintail.quan";
            }else if ((!filter) && (maskfile != "")) {
                tempQuan = tryPath + util.getRootName(util.getSimpleName(baseName)) + "pintail.masked.quan";
            }else if ((filter) && (maskfile != "")) {
                tempQuan = tryPath + util.getRootName(util.getSimpleName(baseName)) + "pintail.filtered." + util.getSimpleName(util.getRootName(fastafile)) + "masked.quan";
            }else if ((filter) && (maskfile == "")) {
                tempQuan = tryPath + util.getRootName(util.getSimpleName(baseName)) + "pintail.filtered." + util.getSimpleName(util.getRootName(fastafile)) + "quan";
            }
            
            ifstream FileTest2(tempQuan.c_str());
            if(FileTest2){
                string line = util.getline(FileTest2);
                bool GoodFile = util.checkReleaseVersion(line, current->getVersion());
                if (GoodFile) {  
                    m->mothurOut("I found " + tempQuan + " in your input file directory. I will use it to save time.\n");  quanfile = tempQuan;  FileTest2.close();	
                }
            }
        }
        
        return 0;
    }
    catch(exception& e) {
        m->errorOut(e, "ChimeraPintailCommand", "lookForShortcutFiles");
        exit(1);
    }
}

/**************************************************************************************************/


